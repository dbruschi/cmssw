import FWCore.ParameterSet.Config as cms
from DQM.TrackingMonitorSource.StandaloneTrackMonitor_cfi import *
from DQM.TrackingMonitorSource.ZEEDetails_cfi import *
from DQM.TrackingMonitorSource.V0Selections_cfi import *
from DQM.TrackingMonitor.V0Monitor_cfi import *

# Primary Vertex Selector
selectedPrimaryVertices = cms.EDFilter("VertexSelector",
                                       src = cms.InputTag('offlinePrimaryVertices'),
                                       # cut = cms.string("!isFake && ndof >= 4 && abs(z) < 24 && abs(position.Rho) < 2.0"),
                                       cut = cms.string(""),
                                       filter = cms.bool(True)
                                       )
# Track Selector
selectedTracks = cms.EDFilter("TrackSelector",
                              src = cms.InputTag('generalTracks'),
                              cut = cms.string("pt > 1.0"),
                              #cut = cms.string(""),
                              filter = cms.bool(True)
                              )

selectedTracksFPixHole = cms.EDFilter("TrackSelector",
                                      src = cms.InputTag('selectedTracks'),
                                      cut = cms.string("eta < -1.4 && phi > 2.5"),
                                      #cut = cms.string(""),
                                      filter = cms.bool(True)
                                      )

# Track Multiplicity Selector
selectedMultiplicityTracks = cms.EDFilter("TrackMultiplicityFilter",
                                          src = cms.InputTag('generalTracks'),
                                          #cut = cms.string("pt > 1.0"),
                                          nmin = cms.untracked.uint32(500),
                                          cut = cms.untracked.string(""),
                                          filter = cms.bool(True)
                                      )

# Track ALCARECO Selection for zerobias
selectedAlcaRecoZBTracks = cms.EDProducer("AlcaRecoTrackSelector",
                                        src = cms.InputTag('generalTracks'),
                                        #cut = cms.string("pt > 0.65 && abs(eta) < 3.5 && p > 1.5 && hitPattern.numberOfAllHits('TRACK_HITS') > 7"),
                                        #cut = cms.string(""),
                                        ptmin = cms.untracked.double(0.65),
                                        pmin = cms.untracked.double(1.5),
                                        etamin = cms.untracked.double(-3.5),
                                        etamax = cms.untracked.double(3.5),
                                        nhits = cms.untracked.uint32(7)
)
'''
# Track ALCARECO Selection for singlemuon
selectedAlcaRecoSMTracks = cms.EDFilter("TrackSelector",
                              src = cms.InputTag('selectedMultiplicityTracks'),
                              cut = cms.string("pt > 2.0 && abs(eta) < 3.5 && p > 1.5 && hitPattern.numberOfAllTrackerHits > 7"),
                              #cut = cms.string(""),
                              filter = cms.bool(True)
                              )
'''
# HLT path selector
hltPathFilter = cms.EDFilter("HLTPathSelector",
                             processName = cms.string("HLT"),
                             verbose = cms.untracked.bool(False),
                             hltPathsOfInterest = cms.vstring("HLT_ZeroBias_v"),
                             triggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                             triggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT")
                             )

# HLT path selector Muon
hltPathFilterMuon = cms.EDFilter("HLTPathSelector",
                             processName = cms.string("HLT"),
                             verbose = cms.untracked.bool(False),
                             hltPathsOfInterest = cms.vstring("HLT_IsoMu24_v"),
                             triggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                             triggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT")
                             )

# HLT path selector Electron
hltPathFilterElectron = cms.EDFilter("HLTPathSelector",
                             processName = cms.string("HLT"),
                             verbose = cms.untracked.bool(False),
                             hltPathsOfInterest = cms.vstring("HLT_Ele32_WPTight_Gsf_v"),
                             triggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                             triggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT")
                             )

# HLT path selector ttbar
hltPathFilterTtbar = cms.EDFilter("HLTPathSelector",
                             processName = cms.string("HLT"),
                             verbose = cms.untracked.bool(False),
                             hltPathsOfInterest = cms.vstring("HLT_Ele32_WPTight_Gsf_v","HLT_IsoMu24_v"),
                             triggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                             triggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT")
                             )

# Z->MuMu event selector
ztoMMEventSelector = cms.EDFilter("ZtoMMEventSelector")
muonTracks = cms.EDProducer("ZtoMMMuonTrackProducer")
muonTracksFPixHole = cms.EDFilter("TrackSelector",
                                  src = cms.InputTag('muonTracks'),
                                  cut = cms.string("eta < -1.4 && phi > 2.5"),
                                  #cut = cms.string(""),
                                  filter = cms.bool(True)
                                  )
# Z->ee event selector
ztoEEEventSelector = cms.EDFilter("ZtoEEEventSelector")
electronTracks = cms.EDProducer("ZtoEEElectronTrackProducer")
electronTracksFPixHole = cms.EDFilter("TrackSelector",
                                      src = cms.InputTag('electronTracks'),
                                      cut = cms.string("eta < -1.4 && phi > 2.5"),
                                      #cut = cms.string(""),
                                      filter = cms.bool(True)
                                      )
#ttbar event selector
ttbarEventSelector = cms.EDFilter("ttbarEventSelector")
ttbarTracks = cms.EDProducer("TtbarTrackProducer")

# Added modules for V0Monitoring
KshortMonitor = v0Monitor.clone()
KshortMonitor.FolderName = "StandaloneTrackMonitor/V0Monitoring/Ks"
KshortMonitor.v0         = "generalV0Candidates:Kshort"
KshortMonitor.histoPSet.massPSet = cms.PSet(nbins = cms.int32 (100),
                                            xmin  = cms.double(0.400),
                                            xmax  = cms.double(0.600))

LambdaMonitor = v0Monitor.clone()
LambdaMonitor.FolderName = "StandaloneTrackMonitor/V0Monitoring/Lambda"
LambdaMonitor.v0 = "generalV0Candidates:Lambda"
LambdaMonitor.histoPSet.massPSet = cms.PSet(nbins = cms.int32(100),
                                            xmin  = cms.double(1.050),
                                            xmax  = cms.double(1.250))
##################
# For MinBias
##################
standaloneTrackMonitorMC = standaloneTrackMonitor.clone(
    puScaleFactorFile = "PileupScaleFactor_316060_wrt_nVertex_ZeroBias.root",
    doPUCorrection    = True,
    isMC              = True
    )
standaloneTrackMonitorFPixHole = standaloneTrackMonitor.clone(
    folderName        = "highPurityTracksFPixHole",
    trackInputTag     = "selectedTracksFPixHole",
    )
standaloneTrackMonitorMCFPixHole = standaloneTrackMonitor.clone(
    folderName        = "highPurityTracksFPixHole",
    trackInputTag     = "selectedTracksFPixHole",
    puScaleFactorFile = "PileupScaleFactor_316060_wrt_nVertex_ZeroBias.root",
    doPUCorrection    = True,
    isMC              = True
    )
standaloneValidationMinbias = cms.Sequence(
    hltPathFilter
    * selectedPrimaryVertices 
#    * selectedMultiplicityTracks  # Use selectedMultiplicityTracks if needed nTracks > desired multiplicity
#    * selectedAlcaRecoZBTracks
    * selectedTracks
    * selectedTracksFPixHole
    * standaloneTrackMonitor
    * standaloneTrackMonitorFPixHole
    * KshortMonitor
    * LambdaMonitor)

standaloneValidationMinbiasMC = cms.Sequence(
    hltPathFilter
    * selectedPrimaryVertices 
#    * selectedMultiplicityTracks  # Use selectedMultiplicityTracks if needed nTracks > desired multiplicity
#    * selectedAlcaRecoZBTracks
    * selectedTracks
    * selectedTracksFPixHole
    * standaloneTrackMonitorMC
    * standaloneTrackMonitorMCFPixHole
    * KshortMonitor
    * LambdaMonitor)

##################
# For V0s in MinBias
##################
standaloneTrackMonitorK0 = standaloneTrackMonitor.clone(
    folderName = "K0Tracks",
    trackInputTag = 'KshortTracks',
    )
KshortTracksFPixHole = cms.EDFilter("TrackSelector",
                                  src = cms.InputTag('KshortTracks'),
                                  cut = cms.string("eta < -1.4 && phi > 2.5"),
                                  #cut = cms.string(""),
                                  filter = cms.bool(True)
                                  )
standaloneTrackMonitorK0FPixHole = standaloneTrackMonitor.clone(
    folderName = "K0TracksFPixHole",
    trackInputTag = 'KshortTracksFPixHole',
    )

standaloneTrackMonitorK0MC = standaloneTrackMonitor.clone(
    folderName = "K0Tracks",
    trackInputTag = 'KshortTracks',
    puScaleFactorFile = "PileupScaleFactor_316082_wrt_nVertex_DYToLL.root",
    doPUCorrection    = True,
    isMC              = True
    )
standaloneTrackMonitorK0MCFPixHole = standaloneTrackMonitor.clone(
    folderName = "K0TracksFPixHole",
    trackInputTag = 'KshortTracksFPixHole',
    puScaleFactorFile = "PileupScaleFactor_316082_wrt_nVertex_DYToLL.root",
    doPUCorrection    = True,
    isMC              = True
    )

standaloneTrackMonitorLambda = standaloneTrackMonitor.clone(
    folderName = "LambdaTracks",
    trackInputTag = 'LambdaTracks',
    )
LambdaTracksFPixHole = cms.EDFilter("TrackSelector",
                                  src = cms.InputTag('LambdaTracks'),
                                  cut = cms.string("eta < -1.4 && phi > 2.5"),
                                  #cut = cms.string(""),
                                  filter = cms.bool(True)
                                  )

standaloneTrackMonitorLambdaFPixHole = standaloneTrackMonitor.clone(
    folderName = "LambdaTracksFPixHole",
    trackInputTag = 'LambdaTracksFPixHole',
    )
standaloneTrackMonitorLambdaMC = standaloneTrackMonitor.clone(
    folderName = "LambdaTracks",
    trackInputTag = 'LambdaTracks',
    puScaleFactorFile = "PileupScaleFactor_316082_wrt_nVertex_DYToLL.root",
    doPUCorrection    = True,
    isMC              = True
    )
standaloneTrackMonitorLambdaMCFPixHole = standaloneTrackMonitor.clone(
    folderName = "LambdaTracksFPixHole",
    trackInputTag = 'LambdaTracksFPixHole',
    puScaleFactorFile = "PileupScaleFactor_316082_wrt_nVertex_DYToLL.root",
    doPUCorrection    = True,
    isMC              = True
    )

standaloneValidationK0s = cms.Sequence(
    hltPathFilter
    * selectedPrimaryVertices
    * KShortEventSelector
    * KshortTracks
    * KshortTracksFPixHole
    * standaloneTrackMonitorK0
    * standaloneTrackMonitorK0FPixHole
    * KshortMonitor)

standaloneValidationK0sMC = cms.Sequence(
    hltPathFilter
    * selectedPrimaryVertices
    * KShortEventSelector
    * KshortTracks
    * KshortTracksFPixHole
    * standaloneTrackMonitorK0MC
    * standaloneTrackMonitorK0MCFPixHole
    * KshortMonitor)

standaloneValidationLambdas = cms.Sequence(
    hltPathFilter
    * selectedPrimaryVertices
    * LambdaEventSelector
    * LambdaTracks
    * LambdaTracksFPixHole
    * standaloneTrackMonitorLambda
    * standaloneTrackMonitorLambdaFPixHole
    * LambdaMonitor)

standaloneValidationLambdasMC = cms.Sequence(
    hltPathFilter
    * selectedPrimaryVertices
    * LambdaEventSelector
    * LambdaTracks
    * LambdaTracksFPixHole
    * standaloneTrackMonitorLambdaMC
    * standaloneTrackMonitorLambdaMCFPixHole
    * LambdaMonitor)

##################
# For ZtoEE
##################
standaloneTrackMonitorElec = standaloneTrackMonitor.clone(
    folderName = "ElectronTracks",
    trackInputTag = 'electronTracks',
    )
standaloneTrackMonitorElecFPixHole = standaloneTrackMonitor.clone(
    folderName = "ElectronTracksFPixHole",
    trackInputTag = 'electronTracksFPixHole',
    )

standaloneTrackMonitorElecMC = standaloneTrackMonitor.clone(
    folderName = "ElectronTracks",
    trackInputTag = 'electronTracks',
    puScaleFactorFile = "PileupScaleFactor_316082_wrt_nVertex_DYToLL.root",
    doPUCorrection    = True,
    isMC              = True
    )
standaloneTrackMonitorElecMCFPixHole = standaloneTrackMonitor.clone(
    folderName = "ElectronTracksFPixHole",
    trackInputTag = 'electronTracksFPixHole',
    puScaleFactorFile = "PileupScaleFactor_316082_wrt_nVertex_DYToLL.root",
    doPUCorrection    = True,
    isMC              = True
    )

ZEEDetailsMC = ZEEDetails.clone(
    puScaleFactorFile = "PileupScaleFactor_316082_wrt_nVertex_DYToLL.root",
    doPUCorrection    = True,
    isMC              = True
    )

standaloneValidationElec = cms.Sequence(
    hltPathFilterElectron
    * selectedTracks
    * selectedTracksFPixHole
    * selectedPrimaryVertices
    * ztoEEEventSelector
    * electronTracks
    * electronTracksFPixHole
    * standaloneTrackMonitorElec
    * standaloneTrackMonitorElecFPixHole   
    * standaloneTrackMonitor
    * standaloneTrackMonitorFPixHole
    * ZEEDetails)
standaloneValidationElecMC = cms.Sequence(
    hltPathFilterElectron
    * selectedTracks
    * selectedTracksFPixHole
    * selectedPrimaryVertices
    * ztoEEEventSelector
    * electronTracks
    * electronTracksFPixHole
    * standaloneTrackMonitorElecMC   
    * standaloneTrackMonitorElecMCFPixHole
    * standaloneTrackMonitorMC
    * standaloneTrackMonitorMCFPixHole
    * ZEEDetailsMC)

##################
# For ZtoMM
##################
standaloneTrackMonitorMuon = standaloneTrackMonitor.clone(
    folderName = "MuonTracks",
    trackInputTag = 'muonTracks',
    )
standaloneTrackMonitorMuonFPixHole = standaloneTrackMonitor.clone(
    folderName = "MuonTracksFPixHole",
    trackInputTag = 'muonTracksFPixHole',
    )
standaloneTrackMonitorMuonMC = standaloneTrackMonitor.clone(
    folderName = "MuonTracks",
    trackInputTag = 'muonTracks',
    puScaleFactorFile = "PileupScaleFactor_316082_wrt_nVertex_DYToLL.root",
    doPUCorrection    = True,
    isMC              = True
    )
standaloneTrackMonitorMuonMCFPixHole = standaloneTrackMonitor.clone(
    folderName = "MuonTracksFPixHole",
    trackInputTag = 'muonTracksFPixHole',
    puScaleFactorFile = "PileupScaleFactor_316082_wrt_nVertex_DYToLL.root",
    doPUCorrection    = True,
    isMC              = True
    )

standaloneValidationMuon = cms.Sequence(
    hltPathFilterMuon
    * selectedTracks
    * selectedTracksFPixHole
    * selectedPrimaryVertices
    * ztoMMEventSelector
    * muonTracks
    * muonTracksFPixHole
    * standaloneTrackMonitorMuon
    * standaloneTrackMonitorMuonFPixHole
    * standaloneTrackMonitor
    * standaloneTrackMonitorFPixHole)

standaloneValidationMuonMC = cms.Sequence(
    hltPathFilterMuon
    * selectedTracks
    * selectedTracksFPixHole
    * selectedPrimaryVertices
    * ztoMMEventSelector
    * muonTracks
    * muonTracksFPixHole
    * standaloneTrackMonitorMuonMC
    * standaloneTrackMonitorMuonMCFPixHole 
    * standaloneTrackMonitorMC
    * standaloneTrackMonitorMCFPixHole)

##################
# For ttbar
##################
standaloneTrackMonitorTTbar = standaloneTrackMonitor.clone(
    folderName = "TTbarTracks",
    trackInputTag = 'ttbarTracks',
    )

standaloneTrackMonitorTTbarMC = standaloneTrackMonitor.clone(
    folderName = "TTbarTracks",
    trackInputTag = 'ttbarTracks',
    puScaleFactorFile = "PileupScaleFactor_316082_wrt_nVertex_DYToLL.root",
    doPUCorrection    = True,
    isMC              = True
    )

standaloneValidationTTbar = cms.Sequence(
    hltPathFilterTtbar
    * selectedPrimaryVertices
    * ttbarEventSelector
    * ttbarTracks
#    * selectedTracks
    * standaloneTrackMonitorTTbar)

standaloneValidationTTbarMC = cms.Sequence(
    hltPathFilterTtbar
    * selectedPrimaryVertices
    * ttbarEventSelector
    * ttbarTracks
#    * selectedTracks
    * standaloneTrackMonitorTTbarMC)
