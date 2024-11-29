import FWCore.ParameterSet.Config as cms


TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
                                                 ComponentName = cms.string('TransientTrackBuilder'),
                                                 MagneticFieldLabel = cms.string("")
                                                 )


