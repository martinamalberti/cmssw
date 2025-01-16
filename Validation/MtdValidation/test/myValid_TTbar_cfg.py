import FWCore.ParameterSet.Config as cms


from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('mtdValidation',Phase2C17I13M9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')
process.load('RecoLocalFastTime.FTLClusterizer.MTDCPEESProducer_cfi')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#Setup FWK for multithreaded
process.options.numberOfThreads = 4
process.options.numberOfStreams = 0
process.options.numberOfConcurrentLuminosityBlocks = 0
process.options.eventSetup.numberOfConcurrentIOVs = 1

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(10),
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:step3.root'
        #Zee
        '/store/relval/CMSSW_14_2_0_pre4/RelValZEE_14/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/998ba580-8d01-40a5-bb10-36821e131033.root',
        '/store/relval/CMSSW_14_2_0_pre4/RelValZEE_14/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/03b9939a-bdac-4d30-905c-2480b29aa970.root',
        '/store/relval/CMSSW_14_2_0_pre4/RelValZEE_14/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/7230a175-7af9-466f-ad67-877d3956f3f2.root',
        '/store/relval/CMSSW_14_2_0_pre4/RelValZEE_14/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/5582d29c-48d0-4c0d-a2d8-d7c468249f18.root',
        '/store/relval/CMSSW_14_2_0_pre4/RelValZEE_14/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/c9d96edb-f25a-44ef-a863-397413699d1a.root',
        '/store/relval/CMSSW_14_2_0_pre4/RelValZEE_14/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/74203ccb-e9e6-4b87-89f8-1d0fd2728574.root',
        '/store/relval/CMSSW_14_2_0_pre4/RelValZEE_14/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/884e8f79-0ae3-4642-aeac-b4c13995e6ae.root',
        '/store/relval/CMSSW_14_2_0_pre4/RelValZEE_14/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/dcc5ef37-f991-40f3-945a-6b56d074b820.root',
        '/store/relval/CMSSW_14_2_0_pre4/RelValZEE_14/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/c358f9b3-7181-49cc-8f74-1dc679083ed5.root',
        '/store/relval/CMSSW_14_2_0_pre4/RelValZEE_14/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/631daaf5-049b-4eee-b2d8-66e286a22fa6.root'
        #ttbar
        #'/store/relval/CMSSW_14_2_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/46671875-91f8-4735-a12a-fa0be47fb1ee.root',
        #'/store/relval/CMSSW_14_2_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/218da3e3-923f-4ae9-b3f6-a3b7d372fc2b.root',
        #'/store/relval/CMSSW_14_2_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/cc049ed3-8720-4f86-b7e1-63f698363ee7.root',
        #'/store/relval/CMSSW_14_2_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/74cebe35-4d74-4be7-b5dd-f405326bd186.root',
        #'/store/relval/CMSSW_14_2_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/41036715-d6d8-471a-8138-57fcdcd1d627.root',
        #'/store/relval/CMSSW_14_2_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/5eaf91c2-8831-4cef-929e-af4b894b7ae0.root',
        #'/store/relval/CMSSW_14_2_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/b2ada0aa-728b-452b-ab24-b6c79b889660.root',
        #'/store/relval/CMSSW_14_2_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/90127bf5-6e0a-49f3-bc44-4ed140145788.root',
        #'/store/relval/CMSSW_14_2_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2590000/7341c2b0-eabd-4a72-bf16-19f88ed9ed03.root'
    )
)

process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)

# --- BTL Validation
process.load("Validation.MtdValidation.btlSimHitsValid_cfi")
process.load("Validation.MtdValidation.btlDigiHitsValid_cfi")
process.load("Validation.MtdValidation.btlLocalRecoValid_cfi")
btlValidation = cms.Sequence(process.btlSimHitsValid + process.btlDigiHitsValid + process.btlLocalRecoValid)

# --- ETL Validation
process.load("Validation.MtdValidation.etlSimHitsValid_cfi")
process.load("Validation.MtdValidation.etlDigiHitsValid_cfi")
process.load("Validation.MtdValidation.etlLocalRecoValid_cfi")
etlValidation = cms.Sequence(process.etlSimHitsValid + process.etlDigiHitsValid + process.etlLocalRecoValid)

# --- Global Validation
process.load("Validation.MtdValidation.mtdTracksValid_cfi")
process.load("Validation.MtdValidation.mtdEleIsoValid_cfi")
process.load("Validation.MtdValidation.vertices4DValid_cfi")

# process.btlDigiHitsValid.optionalPlots = True
# process.etlDigiHitsValid.optionalPlots = True
# process.btlLocalRecoValid.optionalPlots = True
# process.etlLocalRecoValid.optionalPlots = True
# process.mtdTracksValid.optionalPlots = True
# process.vertices4DValid.optionalPlots = True

#process.validation = cms.Sequence(btlValidation + etlValidation + process.mtdTracksValid + process.mtdEleIsoValid + process.vertices4DValid)
process.validation = cms.Sequence(process.mtdEleIsoValid)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.p = cms.Path( process.mix + process.mtdTrackingRecHits + process.validation )
process.endjob_step = cms.EndPath(process.endOfProcess)
process.DQMoutput_step = cms.EndPath( process.DQMoutput )

process.schedule = cms.Schedule( process.p , process.endjob_step , process.DQMoutput_step )
