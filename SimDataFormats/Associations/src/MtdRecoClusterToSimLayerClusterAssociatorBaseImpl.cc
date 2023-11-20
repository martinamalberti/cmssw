#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociatorBaseImpl.h"

namespace reco {
  MtdRecoClusterToSimLayerClusterAssociatorBaseImpl::MtdRecoClusterToSimLayerClusterAssociatorBaseImpl(){};
  MtdRecoClusterToSimLayerClusterAssociatorBaseImpl::~MtdRecoClusterToSimLayerClusterAssociatorBaseImpl(){};

  reco::RecoToSimCollection MtdRecoClusterToSimLayerClusterAssociatorBaseImpl::associateRecoToSim(
      const edm::Handle<FTLClusterCollection> &rCCH, const edm::Handle<MtdSimLayerClusterCollection> &sCCH) const {
    return reco::RecoToSimCollection();
  }

  reco::SimToRecoCollection MtdRecoClusterToSimLayerClusterAssociatorBaseImpl::associateSimToReco(
      const edm::Handle<FTLClusterCollection> &rCCH, const edm::Handle<MtdSimLayerClusterCollection> &sCCH) const {
    return reco::SimToRecoCollection();
  }

}  // namespace reco
