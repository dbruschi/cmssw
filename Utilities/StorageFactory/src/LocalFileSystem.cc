#define _GNU_SOURCE 1
#define _FILE_OFFSET_BITS 64
#include "Utilities/StorageFactory/interface/LocalFileSystem.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <sys/param.h>
#if BSD
#include <sys/statvfs.h>
#include <sys/ucred.h>
#include <sys/mount.h>
#else
#include <mntent.h>
#include <sys/vfs.h>
#endif
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <atomic>
#include <memory>

using namespace edm::storage;

/// Information about file systems on this node.
struct LocalFileSystem::FSInfo {
  char *fsname;               //< file system name
  char *type;                 //< file system type
  char *dir;                  //< mount point directory
  char *origin = nullptr;     //< mount origin
  dev_t dev;                  //< device id
  long fstype;                //< file system id
  double freespc;             //< free space in megabytes
  unsigned local : 1;         //< flag for local device
  unsigned bind : 1;          //< flag for bind mounts
  std::atomic<bool> checked;  //< flag for valid dev, fstype
};

/** Read /proc/filesystems to determine which filesystems are local,
    meaning access latency is tolerably small, and operating system
    buffer cache will likely do a good job at caching file contents
    and accelerate many small file operations reasonably well.

    The /proc list enumerates all filesystems known by the kernel,
    except a few special ones like /dev and /selinux. The ones marked
    as "nodev" have unstable device definition, meaning they are some
    way or another "virtual" file systems.  This labeling is used by
    kernel nfsd to determine which file systems are safe for exporting
    without help (fixing fsid), and turns out to be close enough to
    list of file systems that we can consider to be high-speed local,
    minus a few exceptions.  Everything else we consider "remote" or
    "slow" file systems where application should prefer massive bulk
    streaming I/O for better performance.

    The exceptions to /proc/filesystems list: lustre and fuse file
    systems are forced to remote status. Everything else like NFS,
    AFS, GPFS and various cluster-based systems are already remote. */
int LocalFileSystem::readFSTypes() {
  int ret = 0;

#if __linux__
  constexpr char procfs[] = "/proc/filesystems";
  auto close_ = [](FILE *iFile) { fclose(iFile); };
  std::unique_ptr<FILE, decltype(close_)> fs(fopen(procfs, "r"), close_);
  if (!fs) {
    int nerr = errno;
    edm::LogWarning("LocalFileSystem::readFSTypes()")
        << "Cannot read '" << procfs << "': " << strerror(nerr) << " (error " << nerr << ")";
    return -1;
  }

  ssize_t nread;
  int line = 0;
  auto free_ = [](char **iPtr) { free(*iPtr); };
  while (!feof(fs.get())) {
    char *type = nullptr;
    std::unique_ptr<char *, decltype(free_)> freeType(&type, free_);

    size_t len = 0;
    ++line;

    if ((nread = getdelim(&type, &len, '\t', fs.get())) == -1 && !feof(fs.get())) {
      edm::LogError("LocalFileSystem::readFSTypes()")
          .format("{}:{}: {} ({}; 1)\n", procfs, line, strerror(errno), nread);
      ret = -1;
      break;
    }

    char *fstype = nullptr;
    std::unique_ptr<char *, decltype(free_)> freeFSType(&fstype, free_);
    if ((nread = getdelim(&fstype, &len, '\n', fs.get())) == -1 && !feof(fs.get())) {
      edm::LogError("LocalFileSystem::readFSTypes()")
          .format("{}:{}: {} ({}; 2)\n", procfs, line, strerror(errno), nread);
      ret = -1;
      break;
    }

    if (feof(fs.get())) {
      break;
    }

    if (!strcmp(type, "nodev\t") || !strcmp(fstype, "lustre\n") || !strncmp(fstype, "fuse", 4)) {
      continue;
    }

    assert(nread >= 1);
    fstype[nread - 1] = 0;
    fstypes_.push_back(fstype);
  }
#endif  // __linux__

  return ret;
}

/** Initialise file system description from /etc/mtab info.

    This function saves the information from getmntent(), matching the
    file system type to the known local ones.  It only remembers the
    information from /etc/mtab, so the dev and fstype attributes are
    not yet valid; call statFSInfo() to fill those in.  This avoids
    touching irrelevant filesystems unnecessarily; the file system may
    not be fully functional, or partially offline, or just very slow. */
LocalFileSystem::FSInfo *LocalFileSystem::initFSInfo(void *arg) {
#if BSD
  struct statfs *m = static_cast<struct statfs *>(arg);
  size_t infolen = sizeof(struct FSInfo);
  size_t fslen = strlen(m->f_mntfromname) + 1;
  size_t dirlen = strlen(m->f_mntonname) + 1;
  size_t typelen = strlen(m->f_fstypename) + 1;
  size_t totlen = infolen + fslen + dirlen + typelen;
  FSInfo *i = (FSInfo *)malloc(totlen);
  char *p = (char *)i;
  i->fsname = strncpy(p += infolen, m->f_mntfromname, fslen);
  i->type = strncpy(p += fslen, m->f_fstypename, typelen);
  i->dir = strncpy(p += typelen, m->f_mntonname, dirlen);
  i->dev = m->f_fsid.val[0];
  i->fstype = m->f_type;
  i->freespc = 0;
  i->bind = 0;
  i->origin = nullptr;
  if (m->f_bsize > 0) {
    i->freespc = m->f_bavail;
    i->freespc *= m->f_bsize;
    i->freespc /= 1024. * 1024. * 1024.;
  }
  /* FIXME: This incorrectly says that mounted disk images are local,
     even if it was mounted from a network server. The alternative is
     to walk up the device tree using either a) process IORegistry to
     get the device tree, which lists devices for disk images, and from
     there translate volume uuid to a mount point; b) parse output from
     'hdiutil info -plist' to determine image-path / dev-entry map. */
  i->local = ((m->f_flags & MNT_LOCAL) ? 1 : 0);
  i->checked = 1;
  return i;

#else   // ! BSD
  mntent *m = static_cast<mntent *>(arg);
  size_t infolen = sizeof(struct FSInfo);
  size_t fslen = strlen(m->mnt_fsname) + 1;
  size_t dirlen = strlen(m->mnt_dir) + 1;
  size_t typelen = strlen(m->mnt_type) + 1;
  size_t originlen = strlen(m->mnt_fsname) + 1;
  size_t totlen = infolen + fslen + dirlen + typelen + originlen;
  FSInfo *i = (FSInfo *)malloc(totlen);
  char *p = (char *)i;
  i->fsname = static_cast<char *>(memcpy(p += infolen, m->mnt_fsname, fslen));
  i->type = static_cast<char *>(memcpy(p += fslen, m->mnt_type, typelen));
  i->dir = static_cast<char *>(memcpy(p += typelen, m->mnt_dir, dirlen));
  i->origin = static_cast<char *>(memcpy(p + dirlen, m->mnt_fsname, originlen));
  i->dev = -1;
  i->fstype = -1;
  i->freespc = 0;
  i->local = 0;
  i->checked = false;
  i->bind = strstr(m->mnt_opts, "bind") != nullptr;

  for (size_t j = 0; j < fstypes_.size() && !i->local; ++j)
    if (fstypes_[j] == i->type)
      i->local = 1;
#endif  // BSD

  return i;
}

/** Initialise the list of currently mounted file systems.

    Reads /etc/mtab (or equivalent) to determine all currently mounted
    file systems, and initialises FSInfo structure for them.  It does
    not yet call statFSInfo() on them, so the device and file type ids
    are not yet complete. */
int LocalFileSystem::initFSList() {
#if BSD
  int rc;
  struct statfs *mtab = 0;
  if ((rc = getmntinfo(&mtab, MNT_NOWAIT)) < 0) {
    int nerr = errno;
    edm::LogWarning("LocalFileSystem::initFSList()")
        << "getmntinfo() failed: " << strerror(nerr) << " (error " << nerr << ")";
    return -1;
  }

  fs_.reserve(rc);
  for (int ix = 0; ix < rc; ++ix)
    fs_.push_back(initFSInfo(&mtab[ix]));

  free(mtab);
#else
  const char *const _PATH_MOUNTED_LINUX = "/proc/self/mounts";
  struct mntent *m;
  FILE *mtab = setmntent(_PATH_MOUNTED_LINUX, "r");
  if (!mtab) {
    int nerr = errno;
    edm::LogWarning("LocalFileSystem::initFSList()")
        << "Cannot read '" << _PATH_MOUNTED_LINUX << "': " << strerror(nerr) << " (error " << nerr << ")";
    return -1;
  }

  fs_.reserve(20);
  while ((m = getmntent(mtab)))
    fs_.push_back(initFSInfo(m));

  endmntent(mtab);
#endif

  return 0;
}

/** Figure out file system device and type ids.

    Calls stat() and statfs() on the file system to determine device
    and file system type ids.  These are required to determine if two
    paths are actually on the same file system.

    This function can be called any number of times.  It only does the
    file system check the first time the function is called. */
int LocalFileSystem::statFSInfo(FSInfo *i) const {
  int ret = 0;
  struct stat s;
  struct statfs sfs;

  if (!i->checked) {
    if (lstat(i->dir, &s) < 0) {
      i->checked = true;

      int nerr = errno;
      if (nerr != ENOENT && nerr != EACCES)
        edm::LogWarning("LocalFileSystem::statFSInfo()")
            << "Cannot lstat('" << i->dir << "'): " << strerror(nerr) << " (error " << nerr << ")";
      return -1;
    }

    if (statfs(i->dir, &sfs) < 0) {
      i->checked = true;
      int nerr = errno;
      edm::LogWarning("LocalFileSystem::statFSInfo()")
          << "Cannot statfs('" << i->dir << "'): " << strerror(nerr) << " (error " << nerr << ")";
      return -1;
    }

    i->dev = s.st_dev;
    i->fstype = sfs.f_type;
    if (sfs.f_bsize > 0) {
      i->freespc = sfs.f_bavail;
      i->freespc *= sfs.f_bsize;
      i->freespc /= 1024. * 1024. * 1024.;
    }
    i->checked = true;
  } else if (i->fstype == -1) {
    errno = ENOENT;
    ret = -1;
  }

  return ret;
}

/** Find the file system @a path was mounted from.  The statfs() and
    stat() information for @a path should be in @a sfs and @a s,
    respectively.

    Finds currently mounted file system that @a path is owned by, and
    returns the FSInfo object for it, or null if no matching live file
    system can be found.  If the return value is non-null, then it is
    guaranteed @a path was on that file system.

    A null return value is possible for certain paths which are not on
    any mounted file system (e.g. /dev or /selinux), or if the file
    system is unavailable or some other way dysfunctional, such as
    dead nfs mount or filesystem does not implement statfs().  */
LocalFileSystem::FSInfo *LocalFileSystem::findMount(const char *path,
                                                    struct statfs *sfs,
                                                    struct stat *s,
                                                    std::vector<std::string> &prev_paths) const {
  for (const auto &old_path : prev_paths) {
    if (!strcmp(old_path.c_str(), path)) {
      edm::LogWarning("LocalFileSystem::findMount()") << "Found a loop in bind mounts; stopping evaluation.";
      return nullptr;
    }
  }

  FSInfo *best = nullptr;
  size_t bestlen = 0;
  size_t len = strlen(path);
  for (size_t i = 0; i < fs_.size(); ++i) {
    // First match simply against the file system path.  We don't
    // touch the file system until the path prefix matches.
    // When we have a path prefix match, check the file system if
    //   we don't have a best match candidate yet, OR
    //   this match is longer (more specific) than the previous best OR
    //   this match is the same length and the previous best isn't local
    // The final condition handles cases such as '/' that can appear twice
    // in the file system list, once as 'rootfs' and once as local fs.
    size_t fslen = strlen(fs_[i]->dir);
    if (!strncmp(fs_[i]->dir, path, fslen) &&
        ((fslen == 1 && fs_[i]->dir[0] == '/') || len == fslen || path[fslen] == '/') &&
        (!best || fslen > bestlen || (fslen == bestlen && !best->local))) {
      // Get the file system device and file system ids.
      if (statFSInfo(fs_[i]) < 0)
        return nullptr;

      // Check the path is on the same device / file system.  If this
      // fails, we found a better prefix match on path, but it's the
      // wrong device, so reset our idea of the best match: it can't
      // be the outer mount any more.  Not sure this is the right
      // thing to do with e.g. loop-back or union mounts.
      if (fs_[i]->dev != s->st_dev || fs_[i]->fstype != sfs->f_type) {
        best = nullptr;
        continue;
      }

      // OK this is better than anything else we found so far.
      best = fs_[i];
      bestlen = fslen;
    }
  }
  // In the case of a bind mount, try looking again at the source directory.
  if (best && best->bind && best->origin) {
    struct stat s2;
    struct statfs sfs2;
    std::unique_ptr<char, decltype(std::free) *> fullpath{realpath(best->origin, nullptr), std::free};

    if (!fullpath)
      fullpath.reset(strdup(best->origin));

    if (lstat(fullpath.get(), &s2) < 0) {
      int nerr = errno;
      edm::LogWarning("LocalFileSystem::findMount()") << "Cannot lstat('" << fullpath.get() << "' alias '" << path
                                                      << "'): " << strerror(nerr) << " (error " << nerr << ")";
      return best;
    }

    if (statfs(fullpath.get(), &sfs2) < 0) {
      int nerr = errno;
      edm::LogWarning("LocalFileSystem::findMount()") << "Cannot statfs('" << fullpath.get() << "' alias '" << path
                                                      << "'): " << strerror(nerr) << " (error " << nerr << ")";
      return best;
    }

    prev_paths.push_back(path);
    LocalFileSystem::FSInfo *new_best = findMount(fullpath.get(), &sfs2, &s2, prev_paths);
    return new_best ? new_best : best;
  }

  return best;
}

/** Determine if @a path is on a file system known to be local.

    Returns @c true if the path is definitely known to be local.
    Returns @c false otherwise, including when it's not possible to
    determine anything about the path at all.

    Does not throw exceptions.  If any errors occur, the errors are
    reported as message logger warnings but the actual error is
    swallowed and the function simply returns @c false. */
bool LocalFileSystem::isLocalPath(const std::string &path) const {
  struct stat s;
  struct statfs sfs;
  std::unique_ptr<char, decltype(std::free) *> fullpath{realpath(path.c_str(), nullptr), std::free};

  if (!fullpath)
    fullpath.reset(strdup(path.c_str()));

  if (lstat(fullpath.get(), &s) < 0) {
    int nerr = errno;
    edm::LogWarning("LocalFileSystem::isLocalPath()") << "Cannot lstat('" << fullpath.get() << "' alias '" << path
                                                      << "'): " << strerror(nerr) << " (error " << nerr << ")";
    return false;
  }

  if (statfs(fullpath.get(), &sfs) < 0) {
    int nerr = errno;
    edm::LogWarning("LocalFileSystem::isLocalPath()") << "Cannot statfs('" << fullpath.get() << "' alias '" << path
                                                      << "'): " << strerror(nerr) << " (error " << nerr << ")";
    return false;
  }

  std::vector<std::string> prev_paths;
  FSInfo *m = findMount(fullpath.get(), &sfs, &s, prev_paths);
  return m ? m->local : false;
}

/** Find a writeable directory among @a paths which is known to be
    local and has at least @a minFreeSpace amount of free space in
    gigabytes.

    The @a paths should contain list of relative or absolute candidate
    directories.  If an entry starts with letter "$" then the value of
    that environment variable is used instead; if the value is $TMPDIR
    and the environment variable is empty, "/tmp" is used instead.

    Returns the first path in @a paths which satisfies the criteria,
    expanded to environment variable value if appropriate, resolved
    to full absolute path.  If no suitable path can be found, returns
    an empty string.

    Does not throw exceptions.  If any serious errors occur, the errors
    are reported as message logger warnings but the actual error is
    swallowed and the directory concerned is skipped.  Non-existent
    and inaccessible directories are silently ignored without warning. */
std::pair<std::string, std::string> LocalFileSystem::findCachePath(const std::vector<std::string> &paths,
                                                                   double minFreeSpace) const {
  struct stat s;
  struct statfs sfs;
  std::ostringstream warningst;
  warningst << "Cannot use lazy-download because:\n";

  for (size_t i = 0, e = paths.size(); i < e; ++i) {
    const char *inpath = paths[i].c_str();
    const char *path = inpath;

    if (*path == '$') {
      char *p = std::getenv(path + 1);
      if (p && *p)
        path = p;
      else if (!strcmp(path, "$TMPDIR"))
        path = "/tmp";
    }

    std::unique_ptr<char, decltype(std::free) *> fullpath{realpath(path, nullptr), std::free};
    if (!fullpath)
      fullpath.reset(strdup(path));

#if 0
    std::cerr /* edm::LogInfo("LocalFileSystem") */
      << "Checking if '" << fullpath.get() << "', from '"
      << inpath << "' is valid cache path with "
      << minFreeSpace << " free space" << std::endl;
#endif

    if (lstat(fullpath.get(), &s) < 0) {
      int nerr = errno;
      if (nerr != ENOENT && nerr != EACCES)
        edm::LogWarning("LocalFileSystem::findCachePath()")
            << "Cannot lstat('" << fullpath.get() << "', from '" << inpath << "'): " << strerror(nerr) << " (error "
            << nerr << ")";
      continue;
    }

    if (statfs(fullpath.get(), &sfs) < 0) {
      int nerr = errno;
      edm::LogWarning("LocalFileSystem::findCachePath()")
          << "Cannot statfs('" << fullpath.get() << "', from '" << inpath << "'): " << strerror(nerr) << " (error "
          << nerr << ")";
      continue;
    }

    std::vector<std::string> prev_paths;
    FSInfo *m = findMount(fullpath.get(), &sfs, &s, prev_paths);
#if 0
    std::cerr /* edm::LogInfo("LocalFileSystem") */
      << "Candidate '" << fullpath.get() << "': "
      << "found=" << (m ? 1 : 0)
      << " local=" << (m && m->local)
      << " free=" << (m ? m->freespc : 0)
      << " access=" << access(fullpath.get(), W_OK)
      << std::endl;
#endif

    if (m && m->local && m->freespc >= minFreeSpace && access(fullpath.get(), W_OK) == 0) {
      return std::make_pair(std::string(fullpath.get()), std::string());
    } else if (m) {
      if (!m->local) {
        warningst << "- The mount " << fullpath.get() << " is not local.\n";
      } else if (m->freespc < minFreeSpace) {
        warningst << " - The mount at " << fullpath.get() << " has only " << m->freespc << " GB free; a minumum of "
                  << minFreeSpace << " GB is required.\n";
      } else if (access(fullpath.get(), W_OK)) {
        warningst << " - The process has no permission to write into " << fullpath.get() << "\n";
      }
    }
  }

  std::string warning_str = warningst.str();
  if (!warning_str.empty()) {
    warning_str = warning_str.substr(0, warning_str.size() - 2);
  }

  return std::make_pair(std::string(), std::move(warning_str));
}

/** Initialise local file system status.  */
LocalFileSystem::LocalFileSystem() {
  if (readFSTypes() < 0)
    return;

  if (initFSList() < 0)
    return;
}

/** Free local file system status resources. */
LocalFileSystem::~LocalFileSystem() {
  for (size_t i = 0, e = fs_.size(); i < e; ++i)
    free(fs_[i]);
}
