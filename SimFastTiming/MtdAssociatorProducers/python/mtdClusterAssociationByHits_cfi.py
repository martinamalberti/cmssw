import FWCore.ParameterSet.Config as cms

mtdRecoClusterToSimLayerClusterAssociationByHits = cms.EDProducer('MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer',
   energyCut = cms.double(9999.),
   timeCut   = cms.double(9999.),
   mightGet  = cms.optional.untracked.vstring
)
