import FWCore.ParameterSet.Config as cms

mtdRecoClusterToSimLayerClusterAssociation = cms.EDProducer("MtdRecoClusterToSimLayerClusterAssociatorProducer",
    associator = cms.InputTag('associator'), #?
    mtdSimClustersTag = cms.InputTag('mix','MergedMtdTruthLC'),
    #mtdRecoClustersTag = cms.InputTag('mtdClusters', 'FTLBarrel'),
    btlRecoClustersTag = cms.InputTag('mtdClusters', 'FTLBarrel'),
    etlRecoClustersTag = cms.InputTag('mtdClusters', 'FTLBarrel'),
    btlRecHitsTag = cms.InputTag('mtdRecHits', 'FTLBarrel'),
    etlRecHitsTag = cms.InputTag('mtdRecHits', 'FTLEndcap'),
)
