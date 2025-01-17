import FWCore.ParameterSet.Config as cms

trackProducerFromPatMuons = cms.EDProducer('TrackProducerFromPatMuons',
                                           src = cms.InputTag('slimmedMuons'),
                                           ptMin = cms.double(0.),
                                           innerTrackOnly = cms.bool(True),
                                           mightGet = cms.optional.untracked.vstring
                                           )
