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
    #fileNames = cms.untracked.vstring('file:/eos/cms/store/relval/CMSSW_13_3_0_pre3/RelValZMM_14/GEN-SIM-RECO/131X_mcRun4_realistic_v6_2026D98noPU-v1/2580000/ead207c7-7b5b-4d24-b3d3-762ce9729115.root')
    fileNames = cms.untracked.vstring('file:/eos/cms/store/relval/CMSSW_13_3_0_pre3/RelValZMM_14/GEN-SIM-RECO/PU_131X_mcRun4_realistic_v6_2026D98PU200-v1/2580000/e34950ee-4c9d-43b3-a0e6-c47682e90813.root')
    #fileNames = cms.untracked.vstring('file:/eos/cms/store/relval/CMSSW_13_3_0_pre3/RelValSingleMuPt10/GEN-SIM-RECO/131X_mcRun4_realistic_v6_2026D98noPU-v1/2580000/0df4bc9a-e8da-44dd-b384-36cd32185555.root')
    #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/malberti/MTD/DPG/CMSSW_13_3_0_pre3/mywork/24807.0_SingleMuPt10+2026D98/step3_RAW2DIGI_RECO_RECOSIM_PAT_VALIDATION_DQM.root')
)

process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)

# --- Tracks dumper
process.mtdTracksDumper = cms.EDAnalyzer('MtdTracksDumper',
  inputTagG = cms.InputTag('generalTracks'),
  inputTagT = cms.InputTag('trackExtenderWithMTD'),
  inputTagV = cms.InputTag('offlinePrimaryVertices4D'),
  inputTagH = cms.InputTag('generatorSmeared'),
  SimTag = cms.InputTag('mix', 'MergedTrackTruth'),
  TPtoRecoTrackAssoc = cms.InputTag('trackingParticleRecoTrackAsssociation'),
  tmtd = cms.InputTag('trackExtenderWithMTD', 'generalTracktmtd'),
  sigmatmtd = cms.InputTag('trackExtenderWithMTD', 'generalTracksigmatmtd'),
  t0Src = cms.InputTag('trackExtenderWithMTD', 'generalTrackt0'),
  sigmat0Src = cms.InputTag('trackExtenderWithMTD', 'generalTracksigmat0'),
  trackAssocSrc = cms.InputTag('trackExtenderWithMTD', 'generalTrackassoc'),
  pathLengthSrc = cms.InputTag('trackExtenderWithMTD', 'generalTrackPathLength'),
  btlMatchChi2Tag = cms.InputTag('trackExtenderWithMTD', 'btlMatchChi2'),
  btlMatchTimeChi2Tag = cms.InputTag('trackExtenderWithMTD', 'btlMatchTimeChi2'),
  etlMatchChi2Tag = cms.InputTag('trackExtenderWithMTD', 'etlMatchChi2'),
  etlMatchTimeChi2Tag = cms.InputTag('trackExtenderWithMTD', 'etlMatchTimeChi2'),
  t0SafePID = cms.InputTag('tofPID', 't0safe'),
  sigmat0SafePID = cms.InputTag('tofPID', 'sigmat0safe'),
  sigmat0PID = cms.InputTag('tofPID', 'sigmat0'),
  t0PID = cms.InputTag('tofPID', 't0'),
  trackMVAQual = cms.InputTag('mtdTrackQualityMVA', 'mtdQualMVA'),
  trackMinimumPt = cms.double(0.7),
  trackMaximumBtlEta = cms.double(1.5),
  trackMinimumEtlEta = cms.double(1.6),
  trackMaximumEtlEta = cms.double(3),
  mightGet = cms.optional.untracked.vstring
)


# Output TFile
process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('mtdTracksTree.root'))

process.p = cms.Path(process.mtdTracksDumper)
