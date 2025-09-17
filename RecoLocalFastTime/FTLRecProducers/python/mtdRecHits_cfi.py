import FWCore.ParameterSet.Config as cms

_barrelAlgo = cms.PSet(
    algoName = cms.string("MTDRecHitAlgo"),
    #thresholdToKeep = cms.double(1.),          # MeV
    thresholdToKeep = cms.vdouble(1., 1., 1.5, 1.5, 2.0, 2.0),          # MeV, per RU
    calibrationConstant = cms.double(0.03125), # MeV/pC
)


_endcapAlgo = cms.PSet(
    algoName = cms.string("MTDRecHitAlgo"),
    #thresholdToKeep = cms.double(0.0425),    # MeV
    thresholdToKeep = cms.vdouble(0.0425),    # MeV
    calibrationConstant = cms.double(0.085), # MeV/MIP
)

from Configuration.Eras.Modifier_phase2_etlV4_cff import phase2_etlV4
#phase2_etlV4.toModify(_endcapAlgo, thresholdToKeep = 0.005, calibrationConstant = 0.015 )
phase2_etlV4.toModify(_endcapAlgo, thresholdToKeep = {0.005}, calibrationConstant = 0.015 )

mtdRecHits = cms.EDProducer(
    "MTDRecHitProducer",
    barrel = _barrelAlgo,
    endcap = _endcapAlgo,
    barrelUncalibratedRecHits = cms.InputTag('mtdUncalibratedRecHits:FTLBarrel'),
    endcapUncalibratedRecHits = cms.InputTag('mtdUncalibratedRecHits:FTLEndcap'),
    BarrelHitsName = cms.string('FTLBarrel'),
    EndcapHitsName = cms.string('FTLEndcap'),
)
