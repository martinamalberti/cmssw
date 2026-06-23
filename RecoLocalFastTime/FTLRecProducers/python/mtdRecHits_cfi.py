import FWCore.ParameterSet.Config as cms

_barrelAlgo = cms.PSet(
    algoName = cms.string("MTDRecHitAlgo"),
    thresholdToKeep = cms.double(1.), # [MeV]
    calibrationConstant = cms.double(1.)
)


_endcapAlgo = cms.PSet(
    algoName = cms.string("MTDRecHitAlgo"),
    thresholdToKeep = cms.double(0.0425),    # MeV
    calibrationConstant = cms.double(0.085), # MeV/MIP
)

from Configuration.Eras.Modifier_phase2_etlV4_cff import phase2_etlV4
phase2_etlV4.toModify(_endcapAlgo, thresholdToKeep = 0.005, calibrationConstant = 0.015 )

mtdRecHits = cms.EDProducer(
    "MTDRecHitProducer",
    barrel = _barrelAlgo,
    endcap = _endcapAlgo,
    barrelUncalibratedRecHits = cms.InputTag('mtdUncalibratedRecHits:FTLBarrel'),
    endcapUncalibratedRecHits = cms.InputTag('mtdUncalibratedRecHits:FTLEndcap'),
    BarrelHitsName = cms.string('FTLBarrel'),
    EndcapHitsName = cms.string('FTLEndcap'),
)

from SimFastTiming.FastTimingCommon.mtdDigitizer_cfi import mtdDigitizer
btlRecHitsSoA = cms.EDProducer('btlrechit::BTLRecHitSoAProducer@alpaka',
    baserh = cms.InputTag("btlBaseRecHitsSoA"),
    invLightSpeedLYSO = mtdDigitizer.barrelDigitizer.DeviceSimulation.LightCollectionSlope, # [ns/cm]
    thresholdToKeep = cms.double(0.001), # [GeV]
    calibrationConstant = cms.double(1.), # [GeV / GeV]
    npeSaturationCorr0 = cms.double(-8.54e-06),
    npeSaturationCorr1 = cms.double(1.034),
    npePerGeV = cms.double(2.285e6) # [Npe/GeV]
)

etlRecHitsSoA = cms.EDProducer('etlrechit::ETLRecHitSoAProducer@alpaka',
    baserh = cms.InputTag("etlBaseRecHitsSoA"),
    thresholdToKeep = cms.double(4.25e-5), # [GeV]
    calibrationConstant = cms.double(8.5e-5), # [GeV / ns]
    timeCorr_p0 = cms.double(0.967683), # 0.974683 - 0.007, ad hoc correction for bias from global delay removal
    timeCorr_p1 = cms.double(-0.237274),
    timeCorr_p2 = cms.double(0.021455),
    timeCorr_p3 = cms.double(-0.000727429)
)
