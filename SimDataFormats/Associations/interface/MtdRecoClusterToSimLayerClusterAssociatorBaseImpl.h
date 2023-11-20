#ifndef SimDataFormats_Associations_MtdRecoClusterToSimLayerClusterAssociatorBaseImpl_h
#define SimDataFormats_Associations_MtdRecoClusterToSimLayerClusterAssociatorBaseImpl_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerClusterFwd.h"

namespace reco {

  typedef edm::AssociationMap<
      edm::OneToManyWithQualityGeneric<MtdSimLayerClusterCollection, FTLClusterCollection, std::pair<float, float>>>
      SimToRecoCollection;
  typedef edm::AssociationMap<edm::OneToManyWithQualityGeneric<FTLClusterCollection, MtdSimLayerClusterCollection, float>>
      RecoToSimCollection;
  
  class MtdRecoClusterToSimLayerClusterAssociatorBaseImpl {
  public:
    /// Constructor
    MtdRecoClusterToSimLayerClusterAssociatorBaseImpl();
    /// Destructor
    virtual ~MtdRecoClusterToSimLayerClusterAssociatorBaseImpl();

    /// Associate a MtdRecoCluster to MtdSimLayerClusters
    virtual reco::RecoToSimCollection associateRecoToSim(const edm::Handle<FTLClusterCollection> &rCCH,
							 const edm::Handle<MtdSimLayerClusterCollection> &sCCH) const;
  

    /// Associate a MtdSimLayerClusters to MtdRecoClusters
    virtual reco::SimToRecoCollection associateSimToReco(const edm::Handle<FTLClusterCollection> &rCCH,
							 const edm::Handle<MtdSimLayerClusterCollection> &sCCH) const;
							 
							 
  };
}  // namespace reco

#endif
