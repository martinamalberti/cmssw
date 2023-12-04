import FWCore.ParameterSet.Config as cms


from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('association',Phase2C17I13M9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load("Configuration.Geometry.GeometryExtended2026D95Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')
process.load('RecoLocalFastTime.FTLClusterizer.MTDCPEESProducer_cfi')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

#Setup FWK for multithreaded
process.options.numberOfThreads = 4
process.options.numberOfStreams = 0
process.options.numberOfConcurrentLuminosityBlocks = 0
process.options.eventSetup.numberOfConcurrentIOVs = 1

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(100),
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_13_3_0_pre5/RelValSingleMuPt10/GEN-SIM-RECO/133X_mcRun4_realistic_v1_2026D98noPU-v1/2590000/07ef646d-f099-4df0-8f7e-4331e159f009.root'
    )
)


# -- Association map producer
from SimFastTiming.MtdAssociatorProducers.mtdClusterAssociationByHits_cfi import mtdRecoClusterToSimLayerClusterAssociationByHits as mtdClusterAssociatorByHitsProducer
from SimFastTiming.MtdAssociatorProducers.mtdClusterAssociation_cfi import mtdRecoClusterToSimLayerClusterAssociation

process.mtdClusterAssociatorByHitsProducer = mtdClusterAssociatorByHitsProducer.clone()
process.mtdRecoClusterToSimLayerClusterAssociation  =  mtdRecoClusterToSimLayerClusterAssociation.clone()

process.associationProducers = cms.Task( process.mtdClusterAssociatorByHitsProducer,
                                         process.mtdRecoClusterToSimLayerClusterAssociation)

#process.load('SimFastTiming.MtdAssociatorProducers.mtdClusterAssociationByHits_cfi')
#process.load('SimFastTiming.MtdAssociatorProducers.mtdClusterAssociation_cfi')
#process.associationProducers = cms.Task( process.mtdRecoClusterToSimLayerClusterAssociationByHits,
#                                         process.mtdRecoClusterToSimLayerClusterAssociation)



process.p = cms.Path(process.associationProducers)

process.out = cms.OutputModule("PoolOutputModule", 
        outputCommands = cms.untracked.vstring(
        'keep *_*_*_*',
        'keep mtdRecoClusterToSimLayerClusterAssociation_*_*_*',
        ),
     fileName = cms.untracked.string('OutputWithAsocciationMaps.root')
 )
 

process.ep = cms.EndPath(process.out)
