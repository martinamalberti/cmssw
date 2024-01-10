import FWCore.ParameterSet.Config as cms

mtdSimLayerClusterToTPAssociatorByTrackId = cms.EDProducer('MtdSimLayerClusterToTPAssociatorByTrackIdProducer',
    #btlSimHitsTag = cms.InputTag('mix', 'g4SimHitsFastTimerHitsBarrel'),
    #etlSimHitsTag = cms.InputTag('mix', 'g4SimHitsFastTimerHitsEndcap'),
    mightGet  = cms.optional.untracked.vstring
)
