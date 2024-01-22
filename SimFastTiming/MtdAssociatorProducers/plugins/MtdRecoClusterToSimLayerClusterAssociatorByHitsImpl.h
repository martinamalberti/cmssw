#include <vector>
#include <memory>

#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociator.h"

namespace edm {
  class EDProductGetter;
}

class MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl : public reco::MtdRecoClusterToSimLayerClusterAssociatorBaseImpl {
public:
  explicit MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl(edm::EDProductGetter const &, double, double);

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
  const double energyCut_;
  const double timeCut_;
};

