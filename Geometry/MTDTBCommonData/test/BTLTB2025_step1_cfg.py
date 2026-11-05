# GEN-SIM for the BTL test beam 2025 geometry with active DUT and REF modules (DD4hep),
# adapted from BTLTB2025_step1_DUTbehind_cfg.py (same beam, the module positions are the same as in the 2025 geometry)
#   cmsRun BTLTB2025_step1_cfg.py maxEvents=200 outputFile=BTLTB2025_step1.root
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.maxEvents = 200
options.outputFile = 'BTLTB2025_step1.root'
options.parseArguments()

from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9
from Configuration.ProcessModifiers.dd4hep_cff import dd4hep
process = cms.Process('SIM', Phase2C22I13M9, dd4hep)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Geometry.MTDTBCommonData.GeometryDD4hepExtendedBTLTB2025Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_0T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedBeamProfile_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
process.source = cms.Source("EmptySource")
process.options = cms.untracked.PSet(numberOfThreads = cms.untracked.uint32(1),
                                     numberOfStreams = cms.untracked.uint32(0))

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('generation_step')),
    dataset = cms.untracked.PSet(dataTier = cms.untracked.string('GEN-SIM'), filterName = cms.untracked.string('')),
    fileName = cms.untracked.string('file:' + options.outputFile),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.genstepfilter.triggerConditions = cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_0T', '')

# Same beam as in BTLTB2025_step1_DUTbehind_cfg.py: pions in the transverse plane (eta = 0) at 52 degrees
# with respect to the normal to the DUT surface, crossing first the REF and then the DUT module.
# With eta = 0 and the default Psi the beam profile is in the plane orthogonal to the beam direction:
# BeamMeanY = -91.015 cm gives the vertex at (71.7, -56.0, 0) cm, on the line through the REF and DUT centres.
beamPhi = 0.9075712110370513
process.VtxSmeared.MinEta = 0.
process.VtxSmeared.MaxEta = 0.
process.VtxSmeared.MinPhi = beamPhi
process.VtxSmeared.MaxPhi = beamPhi
process.VtxSmeared.BeamMeanX = 0.
process.VtxSmeared.BeamMeanY = -91.01524204157639
process.VtxSmeared.BeamPosition = 0.
process.VtxSmeared.BeamSigmaX = 2.
process.VtxSmeared.BeamSigmaY = 2.

process.generator = cms.EDProducer("FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(211),
        MinEta = cms.double(0.),
        MaxEta = cms.double(0.),
        MinPhi = cms.double(beamPhi),
        MaxPhi = cms.double(beamPhi),
        MinE   = cms.double(180.00),
        MaxE   = cms.double(180.00)
    ),
    Verbosity       = cms.untracked.int32(0),
    AddAntiParticle = cms.bool(False)
)

# only the MTD sensitive detector is present in the test beam geometry
process.g4SimHits.OnlySDs = ['MtdSensitiveDetector']

process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)
process.schedule = cms.Schedule(process.generation_step, process.genfiltersummary_step, process.simulation_step,
                                process.endjob_step, process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

for path in process.paths:
    getattr(process, path).insert(0, process.generator)

process.MessageLogger.cerr.FwkReport.reportEvery = 50

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
