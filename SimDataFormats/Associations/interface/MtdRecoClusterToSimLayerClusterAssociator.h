#ifndef SimDataFormats_Associations_MtdRecoClusterToSimLayerClusterAssociator_h
#define SimDataFormats_Associations_MtdRecoClusterToSimLayerClusterAssociator_h
// Original Author:  Martina Malberti

// system include files
#include <memory>

// user include files
#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociatorBaseImpl.h"

// forward declarations

namespace reco {

  class MtdRecoClusterToSimLayerClusterAssociator {
  public:
    MtdRecoClusterToSimLayerClusterAssociator(std::unique_ptr<reco::MtdRecoClusterToSimLayerClusterAssociatorBaseImpl>);
    MtdRecoClusterToSimLayerClusterAssociator() = default;
    MtdRecoClusterToSimLayerClusterAssociator(MtdRecoClusterToSimLayerClusterAssociator &&) = default;
    MtdRecoClusterToSimLayerClusterAssociator &operator=(MtdRecoClusterToSimLayerClusterAssociator &&) = default;
    MtdRecoClusterToSimLayerClusterAssociator(const MtdRecoClusterToSimLayerClusterAssociator &) = delete;  // stop default

    ~MtdRecoClusterToSimLayerClusterAssociator() = default;
    const MtdRecoClusterToSimLayerClusterAssociator &operator=(const MtdRecoClusterToSimLayerClusterAssociator &) =
        delete;  // stop default
    
    // ---------- const member functions ---------------------
    /// Associate RecoCluster to MtdSimLayerCluster
    reco::RecoToSimCollectionMtd associateRecoToSim(const edm::Handle<FTLClusterCollection> &rCCH,
						    const edm::Handle<MtdSimLayerClusterCollection> &sCCH,
						    const edm::Handle<FTLRecHitCollection> &btlRecHitsH,
						    const edm::Handle<FTLRecHitCollection> &etlRecHitsH) const {
      return m_impl->associateRecoToSim(rCCH, sCCH, btlRecHitsH, etlRecHitsH);
    };

    /// Associate MtdSimLayerCluster to RecoCluster
    reco::SimToRecoCollectionMtd associateSimToReco(const edm::Handle<FTLClusterCollection> &rCCH,
						    const edm::Handle<MtdSimLayerClusterCollection> &sCCH,
						    const edm::Handle<FTLRecHitCollection> &btlRecHitsH,
                                                    const edm::Handle<FTLRecHitCollection> &etlRecHitsH) const {
      return m_impl->associateSimToReco(rCCH, sCCH, btlRecHitsH, etlRecHitsH);
    };

  private:
    // ---------- member data --------------------------------
    std::unique_ptr<MtdRecoClusterToSimLayerClusterAssociatorBaseImpl> m_impl;
  };
} // namespace reco 

#endif
