import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys
##################################################################
# Print the contents of the RawDataBuffer EDProduct (corresponding to DAQ raw data)
##################################################################
options = VarParsing.VarParsing ('analysis')

process = cms.Process("DUMP")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("file:BtlDigiToRaw.root")
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    #limit = cms.untracked.int32(0)
    limit = cms.untracked.int32(-1)
)
process.MessageLogger.cerr.dumpBTLRawDataBuffer = cms.untracked.PSet(
    #limit = cms.untracked.int32(0)
    limit = cms.untracked.int32(-1)
)

process.maxEvents.input = -1

process.dumpBTLRawDataBuffer = cms.EDAnalyzer("DumpBTLRawDataBuffer",
                                              #minSLinkID = cms.uint32(0),
    #maxSLinkID = cms.uint32(11),
    rawDataBufferTag = cms.InputTag("btlDigiToRaw","","BTLPACKER"),
)

process.path = cms.Path(process.dumpBTLRawDataBuffer)
