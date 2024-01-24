import FWCore.ParameterSet.Config as cms


from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('association',Phase2C17I13M9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load("Configuration.Geometry.GeometryExtended2026D98Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')
process.load('RecoLocalFastTime.FTLClusterizer.MTDCPEESProducer_cfi')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

#Setup FWK for multithreaded
process.options.numberOfThreads = 1
process.options.numberOfStreams = 0
process.options.numberOfConcurrentLuminosityBlocks = 0
process.options.eventSetup.numberOfConcurrentIOVs = 1

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.threshold = 'DEBUG'
process.MessageLogger.debugModules = ['mtdRecoClusterToSimLayerClusterAssociatorByHits','mtdRecoClusterToSimLayerClusterAssociation']
#process.MessageLogger.debugModules = ['*']

process.Timing = cms.Service("Timing",
  #summaryOnly = cms.untracked.bool(False),
  #useJobReport = cms.untracked.bool(True)
)

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(1),
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # SingleMu, noPU        
        #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuPt10/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2590000/1095a786-6fc3-4cd1-a159-89175f5c868a.root'
        #'file:/afs/cern.ch/work/m/malberti/MTD/DPG/CMSSW_14_0_0_pre1/mywork/24807.0_SingleMuPt10+2026D98/step3.root'
        # SingleMu, PU200        
        '/store//relval/CMSSW_14_0_0_pre1/RelValSingleMuPt10/GEN-SIM-RECO/PU_133X_mcRun4_realistic_v1_2026D98PU200-v1/2590000/3c20ae49-f4b0-472f-aa85-3dade4e7a32c.root'
        # SinglePi, noPU        
        #'/store/relval/CMSSW_14_0_0_pre2/RelValSinglePiFlatPt0p7To10/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2590000/25725b7e-b5a1-4c83-ae86-733c4a04c9d0.root'
        # SinglePi, PU200
        #'/store/relval/CMSSW_14_0_0_pre1/RelValSinglePiFlatPt0p7To10/GEN-SIM-RECO/PU_133X_mcRun4_realistic_v1_2026D98PU200-v1/2590000/e81370bc-b36a-4c6b-9427-01a1badfc188.root'        
    )
)

# -- Association maps producers
process.load('SimFastTiming.MtdAssociatorProducers.mtdRecoClusterToSimLayerClusterAssociatorByHits_cfi')
process.load('SimFastTiming.MtdAssociatorProducers.mtdRecoClusterToSimLayerClusterAssociation_cfi')
process.load('SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociatorByTrackId_cfi')
process.load('SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociation_cfi')
associationProducers = cms.Sequence(
    process.mtdRecoClusterToSimLayerClusterAssociatorByHits +
    process.mtdRecoClusterToSimLayerClusterAssociation +
    process.mtdSimLayerClusterToTPAssociatorByTrackId +
    process.mtdSimLayerClusterToTPAssociation
)

process.p = cms.Path(associationProducers)

process.out = cms.OutputModule("PoolOutputModule", 
        outputCommands = cms.untracked.vstring(
        #'keep *_*_*_*',
            'drop *_*_*_*',
            'keep *_*_MergedMtdTruth_*',
            'keep *_*_MergedMtdTruthLC_*',
            'keep TrackingParticles_*_*_*',
            'keep *_mtdClusters_*_*',
            'keep *_mtdRecoClusterToSimLayerClusterAssociation_*_*',
            'keep *_mtdSimLayerClusterToTPAssociation_*_*',
        ),
    #fileName = cms.untracked.string('OutputWithAssociationMaps_SingleMu_noPU.root')
    fileName = cms.untracked.string('OutputWithAssociationMaps_SingleMu_PU200.root')
    #fileName = cms.untracked.string('OutputWithAssociationMaps_SinglePi_noPU.root')
    #fileName = cms.untracked.string('OutputWithAssociationMaps_SinglePi_PU200.root')
)
 

process.ep = cms.EndPath(process.out)
