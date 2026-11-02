import FWCore.ParameterSet.Config as cms

# ETL test beam 2026: three ETL v11 modules along z, standalone world (no BTL)

from Configuration.Geometry.GeometryDD4hep_cff import *
DDDetectorESProducer.confGeomXMLFiles = cms.FileInPath("Geometry/MTDTBCommonData/data/dd4hep/cmsExtendedGeometryETLTB2026.xml")

from Geometry.MTDNumberingBuilder.mtdNumberingGeometry_cff import *
