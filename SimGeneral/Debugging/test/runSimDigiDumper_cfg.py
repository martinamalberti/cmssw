import FWCore.ParameterSet.Config as cms



process = cms.Process("SimDigiDump")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("Configuration.Geometry.GeometryExtended2026D110Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
            #fileNames = cms.untracked.vstring('file:step2.root')
            fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/malberti/MTD/DPG/BTLNumberingScheme/CMSSW_14_2_0_pre3/mywork/29607.0_SingleMuPt10+2026D110/step2.root')
)

process.load("SimGeneral.Debugging.simDigiDumper_cfi")

process.p1 = cms.Path(process.simDigiDumper)
