import FWCore.ParameterSet.Config as cms

mtdSimLayerClusterToTPAssociation = cms.EDProducer("MtdSimLayerClusterToTPAssociatorEDProducer",
    associator = cms.InputTag('mtdSimLayerClusterToTPAssociatorByUniqueId'),
    mtdSimClustersTag = cms.InputTag('mix','MergedMtdTruthLC'),
    trackingParticlesTag = cms.InputTag('mix', 'MergedTrackTruth')
)
