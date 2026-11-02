import FWCore.ParameterSet.Config as cms

# ETL test beam 2026: three ETL v11 modules along z in the CMS skeleton of the MTD test beam simulations,
# with the full BTL v4

from Geometry.MTDTBCommonData.GeometryDD4hepExtendedETLTB2026_cff import *

# timing
from RecoMTD.DetLayers.mtdDetLayerGeometry_cfi import *
from Geometry.MTDGeometryBuilder.mtdParameters_cff import *
from Geometry.MTDNumberingBuilder.mtdNumberingGeometry_cff import *
from Geometry.MTDNumberingBuilder.mtdTopology_cfi import *
from Geometry.MTDGeometryBuilder.mtdGeometry_cfi import *
from Geometry.MTDGeometryBuilder.idealForDigiMTDGeometry_cff import *
mtdGeometry.applyAlignment = False
