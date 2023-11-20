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
    # SingleMuonPt10 - noPU
    #fileNames = cms.untracked.vstring(
    #    'file:/eos/cms/store/relval/CMSSW_13_3_0_pre4/RelValSingleMuPt10/GEN-SIM-RECO/133X_mcRun4_realistic_v1_Standard_13_3_0_pre4-v4/2590000/1985c0d1-ff93-4c75-9451-8235e6aabfef.root',
    #    'file:/eos/cms/store/relval/CMSSW_13_3_0_pre4/RelValSingleMuPt10/GEN-SIM-RECO/133X_mcRun4_realistic_v1_Standard_13_3_0_pre4-v4/2590000/a10dde78-d0b0-45d0-8b14-923d9278272c.root',
    #    'file:/eos/cms/store/relval/CMSSW_13_3_0_pre4/RelValSingleMuPt10/GEN-SIM-RECO/133X_mcRun4_realistic_v1_Standard_13_3_0_pre4-v4/2590000/0f3b7aaf-61bf-47cc-a9f5-971294f2bbcc.root'
                            #)
                            
    # SinglePion - noPU
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/relval/CMSSW_13_3_0_pre4/RelValSinglePiFlatPt0p7To10/GEN-SIM-RECO/133X_mcRun4_realistic_v1_Standard_13_3_0_pre4-v4/2590000/881fab8c-3022-4ba1-919d-ae2e4fdca7e8.root',
        'file:/eos/cms/store/relval/CMSSW_13_3_0_pre4/RelValSinglePiFlatPt0p7To10/GEN-SIM-RECO/133X_mcRun4_realistic_v1_Standard_13_3_0_pre4-v4/2590000/ac536dae-5a0f-4f1f-bd27-3f36a1aa6e11.root',
        'file:/eos/cms/store/relval/CMSSW_13_3_0_pre4/RelValSinglePiFlatPt0p7To10/GEN-SIM-RECO/133X_mcRun4_realistic_v1_Standard_13_3_0_pre4-v4/2590000/06f543fe-4bdf-4535-8fa6-0071a56a8ce9.root'
    )       
)


process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)


# --- Sim Clusters dumper
process.mtdSimClusterDumper = cms.EDAnalyzer('MtdSimClusterDumper',
    btlSimHitsTag = cms.InputTag('mix', 'g4SimHitsFastTimerHitsBarrel'),
    mtdSimLayerClustersTag = cms.InputTag('mix','MergedMtdTruthLC'),
    hitMinEnergy = cms.double(0.),
)


# Output TFile
process.TFileService = cms.Service('TFileService',
    #fileName = cms.string('mtdSimClustersTree_SingleMuonPt10.root')
    fileName = cms.string('mtdSimClustersTree_SinglePiFlat.root')
)

process.p = cms.Path(process.mix + process.mtdSimClusterDumper)
