#include <vector>
#include <memory>

#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToTPAssociator.h"

namespace edm {
  class EDProductGetter;
}

class MtdSimLayerClusterToTPAssociatorByUniqueIdImpl : public reco::MtdSimLayerClusterToTPAssociatorBaseImpl {
public:
  explicit MtdSimLayerClusterToTPAssociatorByUniqueIdImpl(edm::EDProductGetter const &);

  reco::SimToTPCollectionMtd associateSimToTP(
      const edm::Handle<MtdSimLayerClusterCollection> &simClusH,
      const edm::Handle<TrackingParticleCollection>& trackingParticleH) const override;

  reco::TPToSimCollectionMtd associateTPToSim(
      const edm::Handle<MtdSimLayerClusterCollection> &simClusH,
      const edm::Handle<TrackingParticleCollection>& trackingParticleH) const override;
  
private:
  edm::EDProductGetter const *productGetter_;

};

