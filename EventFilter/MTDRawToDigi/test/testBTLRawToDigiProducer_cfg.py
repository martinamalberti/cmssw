import FWCore.ParameterSet.Config as cms

import Configuration.Geometry.defaultPhase2ConditionsEra_cff as _settings
import Geometry.MTDCommonData.defaultMTDConditionsEra_cff as _mtdgeo
_mtdgeo.check_mtdgeo()
_PH2_GLOBAL_TAG, _PH2_ERA = _settings.get_era_and_conditions(_mtdgeo.MTD_DEFAULT_VERSION)
from Configuration.ProcessModifiers.dd4hep_cff import dd4hep

process = cms.Process("BTLUNPACKER",_PH2_ERA, dd4hep)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

process.load('Geometry.MTDCommonData.GeometryDD4hepExtendedRun4MTDDefaultReco_cff')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(3))

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring("file:/afs/cern.ch/work/m/malberti/MTD/DPG/CMSSW_20_0_0_pre1/mywork/34434.0_TTbar_14TeV+Run4D121/step2.root")
                            fileNames = cms.untracked.vstring("file:/afs/cern.ch/work/m/malberti/MTD/DPG/CMSSW_20_0_0_pre1/src/EventFilter/MTDRawToDigi/test/BtlDigiToRaw.root")
)


process.load("Configuration.StandardSequences.Accelerators_cff")

process.load('EventFilter.MTDRawToDigi.BTLReadoutMapESProducer_cfi')
process.load("EventFilter.MTDRawToDigi.BTLElectronicsToDetIdMappingESProducer_cfi")

# EDProducer
process.btlRawToDigiGPU = cms.EDProducer(
    "BTLRawToDigi@alpaka",
    rawDataBufferTag = cms.InputTag("btlDigiToRaw"),    
    # btlDigisLabel    = 'btlDigis'
)


process.out = cms.OutputModule("PoolOutputModule",
                               splitLevel = cms.untracked.int32(0),
                               outputCommands = cms.untracked.vstring('keep *'),
                               fileName = cms.untracked.string('BtlRawToDigi.root')
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.path1 = cms.Path(process.btlRawToDigiGPU)
process.output = cms.EndPath(process.out)
