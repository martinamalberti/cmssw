#include <vector>
#include <map>
#include <unordered_map>
#include <memory>  // shared_ptr

#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociator.h"

namespace edm {
  class EDProductGetter;
}

class MtdRecoClusterToSimLayerClusterAssociatorImpl : public reco::MtdRecoClusterToSimLayerClusterAssociatorImpl {
public:
  MtdRecoClusterToSimLayerClusterAssociatorImpl(edm::EDProductGetter const &);

  reco::RecoToSimCollection associateRecoToSim(
      const edm::Handle<FTLClusterCollection> &rCCH,
      const edm::Handle<MtdSimLayerClusterCollection> &sCCH) const override;

  reco::SimToRecoCollection associateSimToReco(
      const edm::Handle<FTLClusterCollection> &rCCH,
      const edm::Handle<MtdSimLayerClusterCollection> &sCCH) const override;

private:
  edm::EDProductGetter const *productGetter_;
};

