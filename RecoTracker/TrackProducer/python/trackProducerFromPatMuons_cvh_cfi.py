import FWCore.ParameterSet.Config as cms

trackProducerFromPatMuons = cms.EDProducer('TrackProducerFromPatMuons_cvh',
                                           src = cms.InputTag('slimmedMuons'),
                                           ptMin = cms.double(10.),
                                           innerTrackOnly = cms.bool(True),
                                           mightGet = cms.optional.untracked.vstring
                                           )
