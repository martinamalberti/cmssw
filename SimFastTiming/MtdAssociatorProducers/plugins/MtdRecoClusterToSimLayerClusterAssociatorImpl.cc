//
//
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "MtdRecoClusterToSimLayerClusterAssociatorImpl.h"



using namespace reco;
using namespace std;


/* Constructor */

MtdRecoClusterToSimLayerClusterAssociatorImpl::MtdRecoClusterToSimLayerClusterAssociatorImpl(edm::EDProductGetter const& productGetter)
  : productGetter_(&productGetter){}


//
//---member functions
//

reco::RecoToSimCollectionMtd MtdRecoClusterToSimLayerClusterAssociatorImpl::associateRecoToSim(
      const edm::Handle<FTLClusterCollection>& rCCH, const edm::Handle<MtdSimLayerClusterCollection>& sCCH,
      const edm::Handle<FTLRecHitCollection>& btlRecHitsHandle, const edm::Handle<FTLRecHitCollection>& etlRecHitsHandle) const {

  RecoToSimCollectionMtd outputCollection(productGetter_);

  // do stuff to fill the output collection

  // -- get collections
  const auto& recoClusters = *rCCH.product();
  const auto& simClusters  = *sCCH.product();
  const auto& btlRecHits   = *btlRecHitsHandle.product();
  const auto& etlRecHits   = *etlRecHitsHandle.product();

  // -- loop over reco clusters
  for ( const auto detSetClus : *rCCH) {
  
    int simClusIndex = 0;
    float quality = 0;
    
    for (const auto& recoClus : detSetClus) {
      
      BTLDetId clusId = recoClus.id();
      MTDDetId mtdId(clusId);

      // === Clusters in BTL   -- testing only BTL for now.
      if (mtdId.mtdSubDetector() != MTDDetId::BTL)  continue;

      std::vector<uint64_t> recoClusHitIds;

      // -- loop over hits in the reco cluster and find their ids
      for (int ihit = 0; ihit < recoClus.size(); ++ihit) {
	int hit_row = recoClus.minHitRow() + recoClus.hitOffset()[ihit * 2];
	int hit_col = recoClus.minHitCol() + recoClus.hitOffset()[ihit * 2 + 1];
	
	for (auto recHit : btlRecHits) {
	  BTLDetId hitId(recHit.id().rawId());
	  
	  // -- check the hit position
	  if (hitId.mtdSide() != clusId.mtdSide() || hitId.mtdRR() != clusId.mtdRR() || recHit.row() != hit_row || recHit.column() != hit_col)
	    continue;
	  
	  // -- check the hit energy and time
	  if (recHit.energy() != recoClus.hitENERGY()[ihit] || recHit.time() != recoClus.hitTIME()[ihit])
	    continue;
	  
	  recoClusHitIds.push_back(hitId);
	}
      } // end loop over hits in reco cluster

      // -- loop over sim clusters and if this reco clus shares some hits
      for (const auto& simClus  : simClusters){
	simClusIndex++;
	std::vector<std::pair<uint64_t, float>> hitsAndFrac = simClus.hits_and_fractions();
	std::vector<uint64_t> simClusHitIds(hitsAndFrac.size());
	std::transform(hitsAndFrac.begin(), hitsAndFrac.end(), simClusHitIds.begin(), [](const std::pair<int, float>& pair) {
											return pair.first;});
	std::vector<uint64_t> sharedHitIds;
	std::set_intersection(recoClusHitIds.begin(), recoClusHitIds.end(), simClusHitIds.begin(), simClusHitIds.end(), std::back_inserter(sharedHitIds));
	if (sharedHitIds.empty()) continue;
	quality = sharedHitIds.size()/recoClusHitIds.size();
      }

      auto recoClusterRef = makeRefTo(rCCH, &recoClus);

      // -- if they share some hits fill the output collection
      //if ( simClusIndex != -1){
	//edm::Ref<FTLClusterCollection> recoClusterRef = edm::Ref<FTLClusterCollection>(rCCH, recoClus);
	//edm::Ref<MtdSimLayerClusterCollection> simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(sCCH, simClusters[simClusIndex]);
	//auto simClusterRef = simClusters[simClusIndex];
	//outputCollection.insert(recoClusterRef, std::make_pair(simClusterRef, quality));
	//outputCollection.insert(recoClusterRef, simClusterRef, quality);  
      //}

	  
    }// -- end loop over reco clus

    /// BOH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( simClusIndex != -1){
      auto recoClusterRef = makeRefTo(rCCH, detSetClus);
      //edm::Ref<FTLClusterCollection> recoClusterRef = edm::Ref<FTLClusterCollection>(recoClusters, DetSetClus);
      //edm::Ref<MtdSimLayerClusterCollection> simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusters, simClusIndex);
      // -- insert vuole degli edm::Ref???????????
      outputCollection.insert(recoClusterRef, simClusterRef, quality);   
    }

    
  }// -- end loop over detsetclus
  
  return outputCollection;

}



reco::SimToRecoCollectionMtd MtdRecoClusterToSimLayerClusterAssociatorImpl::associateSimToReco(
      const edm::Handle<FTLClusterCollection>& rCCH, const edm::Handle<MtdSimLayerClusterCollection>& sCCH,
      const edm::Handle<FTLRecHitCollection>& btlRecHitsHandle, const edm::Handle<FTLRecHitCollection>& etlRecHitsHandle) const {

  SimToRecoCollectionMtd outputCollection(productGetter_);

  // do stuff to fill the output collection
  


  return outputCollection;
}
