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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

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
        '/store/relval/CMSSW_14_0_0_pre1/RelValSingleMuPt10/GEN-SIM-RECO/133X_mcRun4_realistic_v1_2026D98noPU-v1/2590000/1a275880-2806-45d6-af7a-403f6d6fc19b.root'
        #'/store/relval/CMSSW_14_0_0_pre1/RelValSinglePiFlatPt0p7To10/GEN-SIM-RECO/133X_mcRun4_realistic_v1_2026D98noPU-v1/2590000/475f5113-6436-4ee4-b266-c6ac0b527126.root'
    )
)

process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)

# -- Association map producer
process.load('SimFastTiming.MtdAssociatorProducers.mtdRecoClusterToSimLayerClusterAssociatorByHits_cfi')
process.load('SimFastTiming.MtdAssociatorProducers.mtdRecoClusterToSimLayerClusterAssociation_cfi')
process.load('SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociatorByTrackId_cfi')
process.load('SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociation_cfi')
associationProducers = cms.Sequence( process.mtdRecoClusterToSimLayerClusterAssociatorByHits + process.mtdRecoClusterToSimLayerClusterAssociation + process.mtdSimLayerClusterToTPAssociatorByTrackId + process.mtdSimLayerClusterToTPAssociation)


process.p = cms.Path(process.mix + associationProducers)

process.out = cms.OutputModule("PoolOutputModule", 
        outputCommands = cms.untracked.vstring(
        'keep *_*_*_*',
        'keep mtdRecoClusterToSimLayerClusterAssociation_*_*_*',
        'keep mtdSimLayerClusterToTPAssociation_*_*_*',
        ),
    fileName = cms.untracked.string('OutputWithAssociationMaps.root')
    #fileName = cms.untracked.string('OutputWithAssocciationMaps_SinglePi.root')
)
 

process.ep = cms.EndPath(process.out)
