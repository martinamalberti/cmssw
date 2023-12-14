#ifndef SimDataFormats_Associations_MtdSimLayerClusterToRecoClusterAssociationMap_h
#define SimDataFormats_Associations_MtdSimLayerClusterToRecoClusterAssociationMap_h

#include "DataFormats/Provenance/interface/ProductID.h"
#include "DataFormats/Common/interface/HandleBase.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerClusterFwd.h"


#include <vector>
#include <utility>
#include <algorithm>

/**
 * Maps MtdSimLayerCluserRef to FTLClusterRef
 *
 */
class MtdSimLayerClusterToRecoClusterAssociationMap {
public:

  using key_type = MtdSimLayerClusterRef;
  using mapped_type = FTLClusterRef;
  using value_type = std::pair<key_type, mapped_type>;
  using map_type = std::vector<value_type>;
  
  /// Constructor
  MtdSimLayerClusterToRecoClusterAssociationMap();
  /// Destructor
  ~MtdSimLayerClusterToRecoClusterAssociationMap();
  
  void emplace_back(const MtdSimLayerClusterRef& simClus, const FTLClusterRef& recoClus) {
    map_.emplace_back(simClus, recoClus);
  }


private:
  map_type map_;

};



#endif
