import FWCore.ParameterSet.Config as cms

mtdRecoClusterToSimLayerClusterAssociation = cms.EDProducer("MtdRecoClusterToSimLayerClusterAssociatorEDProducer",
    associator = cms.InputTag('mtdClusterAssociatorByHitsProducer'),
    mtdSimClustersTag = cms.InputTag('mix','MergedMtdTruthLC'),
    btlRecoClustersTag = cms.InputTag('mtdClusters', 'FTLBarrel'),
    etlRecoClustersTag = cms.InputTag('mtdClusters', 'FTLEndcap'),
    btlRecHitsTag = cms.InputTag('mtdRecHits', 'FTLBarrel'),
    etlRecHitsTag = cms.InputTag('mtdRecHits', 'FTLEndcap'),
)
