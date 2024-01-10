#include <vector>
#include <memory>

#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeomUtil.h"

#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToTPAssociator.h"

namespace edm {
  class EDProductGetter;
}

class MtdSimLayerClusterToTPAssociatorByTrackIdImpl : public reco::MtdSimLayerClusterToTPAssociatorBaseImpl {
public:
  explicit MtdSimLayerClusterToTPAssociatorByTrackIdImpl(edm::EDProductGetter const &
							 //,
							 //const edm::Handle<CrossingFrame<PSimHit>> &,
							 //const edm::Handle<CrossingFrame<PSimHit>> &,
							 //mtd::MTDGeomUtil &
							 );
  
  reco::SimToTPCollectionMtd associateSimToTP(
      const edm::Handle<MtdSimLayerClusterCollection> &simClusH,
      const edm::Handle<TrackingParticleCollection>& trackingParticleH) const override;

  reco::TPToSimCollectionMtd associateTPToSim(
      const edm::Handle<MtdSimLayerClusterCollection> &simClusH,
      const edm::Handle<TrackingParticleCollection>& trackingParticleH) const override;
  
private:

  edm::EDProductGetter const *productGetter_;
  //edm::Handle<CrossingFrame<PSimHit>> btlSimHitsH_;
  //edm::Handle<CrossingFrame<PSimHit>> etlSimHitsH_;
  //mtd::MTDGeomUtil geomTools_;

};

