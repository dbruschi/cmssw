import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
process = cms.Process('RECO2',Run2_2016)

options = VarParsing.VarParsing('analysis')
options.register('globalTag',
                 "auto:run2_mc", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string, # string, int, or float
                 "input file name")
options.parseArguments()

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
    
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles),
                            secondaryFileNames = cms.untracked.vstring()
                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.options = cms.untracked.PSet()

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('step1_RECO.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition
process.RECOSIMoutput.outputCommands = cms.untracked.vstring("keep *_myRefittedTracks_*_*")

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, '')

# track refit stuff
from TrackPropagation.Geant4e.geantRefit_cff import geopro
process.load("TrackPropagation.Geant4e.geantRefit_cff")
from RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff import *
from RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff import *

import RecoTracker.TrackProducer.trackProducerFromPatMuons_cvh_cfi
process.tracksFromMuons  = RecoTracker.TrackProducer.trackProducerFromPatMuons_cvh_cfi.trackProducerFromPatMuons.clone(
    src = "slimmedMuons",
    innerTrackOnly = True
)

process.trackrefit = cms.EDProducer('ResidualGlobalCorrectionMakerG4e',
                                    src = cms.InputTag("tracksFromMuons"),
                                    fitFromGenParms = cms.bool(False),
                                    fitFromSimParms = cms.bool(False),
                                    fillTrackTree = cms.bool(True),
                                    fillGrads = cms.bool(False),
                                    fillJac = cms.bool(False),
                                    fillRunTree = cms.bool(False),
                                    doGen = cms.bool(False),
                                    doSim = cms.bool(False),
                                    requireGen = cms.bool(False),
                                    doMuons = cms.bool(False),
                                    doMuonAssoc = cms.bool(True),
                                    doTrigger = cms.bool(False),
                                    doRes = cms.bool(False),
                                    bsConstraint = cms.bool(False),
                                    applyHitQuality = cms.bool(True),
                                    useIdealGeometry = cms.bool(False),
                                    corFiles = cms.vstring(),
                                    MagneticFieldLabel = cms.string(""),
                                    )

process.trackrefitideal = cms.EDProducer('ResidualGlobalCorrectionMakerG4e',
                                         src = cms.InputTag("tracksFromMuons"),
                                         fitFromGenParms = cms.bool(False),
                                         fitFromSimParms = cms.bool(False),
                                         fillTrackTree = cms.bool(True),
                                         fillGrads = cms.bool(False),
                                         fillJac = cms.bool(False),
                                         fillRunTree = cms.bool(False),
                                         doGen = cms.bool(False),
                                         doSim = cms.bool(False),
                                         requireGen = cms.bool(False),
                                         doMuons = cms.bool(False),
                                         doMuonAssoc = cms.bool(True),
                                         doTrigger = cms.bool(False),
                                         doRes = cms.bool(False),
                                         bsConstraint = cms.bool(False),
                                         applyHitQuality = cms.bool(True),
                                         useIdealGeometry = cms.bool(True),
                                         corFiles = cms.vstring(),
                                         MagneticFieldLabel = cms.string(""),
                                         )

process.mergedGlobalIdxs = cms.EDProducer("GlobalIdxProducer",
                                          src0 = cms.InputTag("trackrefit", "globalIdxs"),
                                          src1 = cms.InputTag("trackrefitideal", "globalIdxs")
                                          )

# import RecoTracker.TrackProducer.TrackRefitter_cfi
# process.myRefittedTracks = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone(
#     src = 'tracksFromMuons',
#     NavigationSchool = '',
#     Fitter = 'FlexibleKFFittingSmoother'
# )
                                         
# Path and EndPath definitions
process.reconstruction_step = cms.Path(process.geopro+process.tracksFromMuons*process.trackrefitideal, #process.myRefittedTracks, 
                                       cms.Task(process.TTRHBuilderAngleAndTemplate, process.templates, process.SiPixelTemplateStoreESProducer))
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step)
