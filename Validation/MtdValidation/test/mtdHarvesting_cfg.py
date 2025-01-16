import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('mtdHarvesting',Phase2C17I13M9)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(-1),
)

# Input source
process.source = cms.Source("DQMRootSource",
                            #fileNames = cms.untracked.vstring('file:step3_inDQM.root')
fileNames = cms.untracked.vstring(#'file:/eos/cms/store/user/malberti/RelValTTbar_14TeV/test_mtdEleIsoValidation_TTbar_PU/241220_170810/0000/step3_inDQM_2.root',
                                  #'file:/eos/cms/store/user/malberti/RelValTTbar_14TeV/test_mtdEleIsoValidation_TTbar_PU/241220_170810/0000/step3_inDQM_5.root',
                                  #'file:/eos/cms/store/user/malberti/RelValTTbar_14TeV/test_mtdEleIsoValidation_TTbar_PU/241220_170810/0000/step3_inDQM_6.root',
                                  #'file:/eos/cms/store/user/malberti/RelValTTbar_14TeV/test_mtdEleIsoValidation_TTbar_PU/241220_170810/0000/step3_inDQM_7.root',
                                  #'file:/eos/cms/store/user/malberti/RelValTTbar_14TeV/test_mtdEleIsoValidation_TTbar_PU/241220_170810/0000/step3_inDQM_8.root',
                                  #'file:/eos/cms/store/user/malberti/RelValTTbar_14TeV/test_mtdEleIsoValidation_TTbar_PU/241220_170810/0000/step3_inDQM_9.root'

                                  'file:/eos/cms/store/user/malberti/RelValZEE_14/test_mtdEleIsoValidation_Zee_PU/241220_170728/0000/step3_inDQM_3.root',
                                  'file:/eos/cms/store/user/malberti/RelValZEE_14/test_mtdEleIsoValidation_Zee_PU/241220_170728/0000/step3_inDQM_5.root',
                                  'file:/eos/cms/store/user/malberti/RelValZEE_14/test_mtdEleIsoValidation_Zee_PU/241220_170728/0000/step3_inDQM_8.root',
                                  'file:/eos/cms/store/user/malberti/RelValZEE_14/test_mtdEleIsoValidation_Zee_PU/241220_170728/0000/step3_inDQM_9.root'
)
)

# Path and EndPath definitions

process.edmtome_step = cms.Path(process.EDMtoME)
process.dqmsave_step = cms.Path(process.DQMSaver)

# --- PostProcessing

process.load("Validation.MtdValidation.btlSimHitsPostProcessor_cfi")
process.load("Validation.MtdValidation.btlLocalRecoPostProcessor_cfi")
process.load("Validation.MtdValidation.MtdTracksPostProcessor_cfi")
process.load("Validation.MtdValidation.MtdEleIsoPostProcessor_cfi")
process.load("Validation.MtdValidation.Primary4DVertexPostProcessor_cfi")

process.harvesting = cms.Sequence(process.btlSimHitsPostProcessor + process.btlLocalRecoPostProcessor + process.MtdTracksPostProcessor + process.MtdEleIsoPostProcessor + process.Primary4DVertexPostProcessor)

process.p = cms.Path( process.harvesting )

process.schedule = cms.Schedule( process.edmtome_step , process.p , process.dqmsave_step )
