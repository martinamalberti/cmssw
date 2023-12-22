import FWCore.ParameterSet.Config as cms

mtdRecoClusterToSimLayerClusterAssociatorByHits = cms.EDProducer('MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer',
    btlRecHitsTag = cms.InputTag('mtdRecHits', 'FTLBarrel'),
    etlRecHitsTag = cms.InputTag('mtdRecHits', 'FTLEndcap'),
    energyCut = cms.double(5.),
    timeCut   = cms.double(10.),
    mightGet  = cms.optional.untracked.vstring
)
