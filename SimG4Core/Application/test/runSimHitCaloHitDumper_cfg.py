import FWCore.ParameterSet.Config as cms

process = cms.Process("SimCaloHitDump")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/malberti/MTD/DPG/BTLNumberingScheme/CMSSW_14_2_0_pre3/mywork/29607.0_SingleMuPt10+2026D110/step1.root')
)

process.load("SimG4Core.Application.simHitCaloHitDumper_cfi")

process.p1 = cms.Path(process.simHitCaloHitDumper)
