import FWCore.ParameterSet.Config as cms

# readout map (DetId -> electronics ids) needed by the packer
from EventFilter.MTDRawToDigi.BTLReadoutMapESProducer_cfi import BTLReadoutMapESProducer

# BTL packer: BTL digis -> RawDataBuffer (btlDigiToRaw_cfi is generated from the fillDescriptions of BTLDigiToRaw)
from EventFilter.MTDRawToDigi.btlDigiToRaw_cfi import btlDigiToRaw

# TODO: add premixing (premix_stage2: BTL digis from mixData), once the MTDBarrel digis are supported in the premixing
