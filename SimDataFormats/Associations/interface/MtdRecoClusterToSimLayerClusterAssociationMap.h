#ifndef SimDataFormats_Associations_MtdRecoClusterToSimLayerClusterAssociationMap_h
#define SimDataFormats_Associations_MtdRecoClusterToSimLayerClusterAssociationMap_h

#include "DataFormats/Provenance/interface/ProductID.h"
#include "DataFormats/Common/interface/HandleBase.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerClusterFwd.h"


#include <vector>
#include <utility>
#include <algorithm>

/**
 * Maps FTLCluserRef to SimLayerClusterRef
 *
 */
class MtdRecoClusterToSimLayerClusterAssociationMap {
public:

  using key_type = FTLClusterRef;
  using mapped_type = MtdSimLayerClusterRef;
  using value_type = std::pair<key_type, mapped_type>;
  using map_type = std::vector<value_type>;
  
  /// Constructor
  MtdRecoClusterToSimLayerClusterAssociationMap();
  /// Destructor
  ~MtdRecoClusterToSimLayerClusterAssociationMap();
  
  void emplace_back(const FTLClusterRef& recoClus, const MtdSimLayerClusterRef& simClus) {
    map_.emplace_back(recoClus, simClus);
  }


private:
  map_type map_;

};



#endif
