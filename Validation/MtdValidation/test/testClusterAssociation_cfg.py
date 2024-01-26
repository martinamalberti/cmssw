import FWCore.ParameterSet.Config as cms


from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('mtdValidation',Phase2C17I13M9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load("Configuration.Geometry.GeometryExtended2026D98Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')
process.load('RecoLocalFastTime.FTLClusterizer.MTDCPEESProducer_cfi')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#Setup FWK for multithreaded
process.options.numberOfThreads = 4
process.options.numberOfStreams = 0
process.options.numberOfConcurrentLuminosityBlocks = 0
process.options.eventSetup.numberOfConcurrentIOVs = 1

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(100),
)

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/malberti/MTD/DPG/CMSSW_14_0_0_pre2/src/SimFastTiming/MtdAssociatorProducers/test/OutputWithAssociationMaps_SingleMu_noPU.root')
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/malberti/MTD/DPG/CMSSW_14_0_0_pre2/src/SimFastTiming/MtdAssociatorProducers/test/OutputWithAssociationMaps_SingleMu_PU200.root')
    #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/malberti/MTD/DPG/CMSSW_14_0_0_pre2/src/SimFastTiming/MtdAssociatorProducers/test/OutputWithAssociationMaps_SinglePi_noPU.root')
    #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/malberti/MTD/DPG/CMSSW_14_0_0_pre2/src/SimFastTiming/MtdAssociatorProducers/test/OutputWithAssociationMaps_SinglePi_PU200.root')

    #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/malberti/MTD/DPG/CMSSW_14_0_0_pre2/src/SimFastTiming/MtdAssociatorProducers/test/OutputWithAssociationMaps_SingleMu_noPU_new.root')
    #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/malberti/MTD/DPG/CMSSW_14_0_0_pre2/src/SimFastTiming/MtdAssociatorProducers/test/OutputWithAssociationMaps_SinglePi_noPU_new.root')

                            )


process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)


# --- Sim Clusters dumper
process.clusterAssociation = cms.EDAnalyzer('TestClusterAssociation',
    btlRecoClustersTag = cms.InputTag('mtdClusters', 'FTLBarrel'),
    etlRecoClustersTag = cms.InputTag('mtdClusters', 'FTLEndcap'),
    r2sAssociationMapTag = cms.InputTag('mtdRecoClusterToSimLayerClusterAssociation'),
    s2rAssociationMapTag = cms.InputTag('mtdRecoClusterToSimLayerClusterAssociation'),
    sim2TPAssociationMapTag = cms.InputTag('mtdSimLayerClusterToTPAssociation'),
    tp2SimAssociationMapTag = cms.InputTag('mtdSimLayerClusterToTPAssociation'),
    mtdSimLayerClustersTag = cms.InputTag('mix','MergedMtdTruthLC'),
    trkHitTag = cms.InputTag('mtdTrackingRecHits')
)


# Output TFile
process.TFileService = cms.Service('TFileService',
   #fileName = cms.string('clusterAssociationTest_SingleMu_noPU.root')
   fileName = cms.string('clusterAssociationTest_SingleMu_PU200.root')
   #fileName = cms.string('clusterAssociationTest_SinglePi_noPU.root')
   #fileName = cms.string('clusterAssociationTest_SinglePi_PU200.root')

   #fileName = cms.string('clusterAssociationTest_SingleMu_noPU_new.root')
   #fileName = cms.string('clusterAssociationTest_SinglePi_noPU_new.root')

                                   )

process.p = cms.Path(process.mix + process.clusterAssociation)
