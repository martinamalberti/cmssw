#include <vector>
#include <map>
#include <unordered_map>
#include <memory>  // shared_ptr

#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociator.h"

namespace edm {
  class EDProductGetter;
}

class MtdRecoClusterToSimLayerClusterAssociatorImpl : public reco::MtdRecoClusterToSimLayerClusterAssociatorBaseImpl {
public:
  MtdRecoClusterToSimLayerClusterAssociatorImpl(edm::EDProductGetter const &);

  reco::RecoToSimCollectionMtd associateRecoToSim(
      const edm::Handle<FTLClusterCollection> &rCCH,
      const edm::Handle<MtdSimLayerClusterCollection> &sCCH,
      const edm::Handle<FTLRecHitCollection>& btlRecHitsHandle,
      const edm::Handle<FTLRecHitCollection>& etlRecHitsHandle) const override;

  reco::SimToRecoCollectionMtd associateSimToReco(
      const edm::Handle<FTLClusterCollection> &rCCH,
      const edm::Handle<MtdSimLayerClusterCollection> &sCCH,
      const edm::Handle<FTLRecHitCollection>& btlRecHitsHandle,
      const edm::Handle<FTLRecHitCollection>& etlRecHitsHandle) const override;

private:
  edm::EDProductGetter const *productGetter_;
};

