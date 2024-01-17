#include <vector>
#include <memory>  // shared_ptr

#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeomUtil.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociator.h"

namespace edm {
  class EDProductGetter;
}

class MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl : public reco::MtdRecoClusterToSimLayerClusterAssociatorBaseImpl {
public:
  explicit MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl(edm::EDProductGetter const &,
							       //const edm::Handle<FTLRecHitCollection>&,
							       //const edm::Handle<FTLRecHitCollection>&,
							       //const MTDTopology*,
							       double,
							       double);

  reco::RecoToSimCollectionMtd associateRecoToSim(
      const edm::Handle<FTLClusterCollection> &btlRecoClusH,
      const edm::Handle<FTLClusterCollection> &etlRecoClusH,
      const edm::Handle<MtdSimLayerClusterCollection> &simClusH) const override;

  reco::SimToRecoCollectionMtd associateSimToReco(
      const edm::Handle<FTLClusterCollection> &btlRecoClusH,
      const edm::Handle<FTLClusterCollection> &etlRecoClusH,
      const edm::Handle<MtdSimLayerClusterCollection> &simClusH) const override;
  
private:
  edm::EDProductGetter const *productGetter_;
  //const edm::Handle<FTLRecHitCollection> btlRecHitsH_;
  //const edm::Handle<FTLRecHitCollection> etlRecHitsH_;
  //const MTDTopology* topo_;
  const double energyCut_;
  const double timeCut_;
};

