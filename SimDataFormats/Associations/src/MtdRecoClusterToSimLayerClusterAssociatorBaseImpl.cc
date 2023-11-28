#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociatorBaseImpl.h"

namespace reco {
  MtdRecoClusterToSimLayerClusterAssociatorBaseImpl::MtdRecoClusterToSimLayerClusterAssociatorBaseImpl(){};
  MtdRecoClusterToSimLayerClusterAssociatorBaseImpl::~MtdRecoClusterToSimLayerClusterAssociatorBaseImpl(){};

  reco::RecoToSimCollectionMtd MtdRecoClusterToSimLayerClusterAssociatorBaseImpl::associateRecoToSim(
      const edm::Handle<FTLClusterCollection> &btlRecoClusH, const edm::Handle<FTLClusterCollection> &etlRecoClusH,
      const edm::Handle<MtdSimLayerClusterCollection> &simClusH,
      const edm::Handle<FTLRecHitCollection> &btlRecHitsH, const edm::Handle<FTLRecHitCollection> &etlRecHitsH) const {
    return reco::RecoToSimCollectionMtd();
  }

  reco::SimToRecoCollectionMtd MtdRecoClusterToSimLayerClusterAssociatorBaseImpl::associateSimToReco(
      const edm::Handle<FTLClusterCollection> &btlRecoClusH, const edm::Handle<FTLClusterCollection> &etlRecoClusH,
      const edm::Handle<MtdSimLayerClusterCollection> &simClusH,
      const edm::Handle<FTLRecHitCollection> &btlRecHitsH, const edm::Handle<FTLRecHitCollection> &etlRecHitsH) const {
    return reco::SimToRecoCollectionMtd();
  }

}  // namespace reco
