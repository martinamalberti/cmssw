//
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "MtdSimLayerClusterToTPAssociatorByTrackIdImpl.h"

#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

#include "SimFastTiming/MtdAssociatorProducers/interface/MtdAssociatorTools.h"

using namespace reco;
using namespace std;


/* Constructor */

MtdSimLayerClusterToTPAssociatorByTrackIdImpl::MtdSimLayerClusterToTPAssociatorByTrackIdImpl(edm::EDProductGetter const& productGetter,
											       const edm::Handle<CrossingFrame<PSimHit>> &btlSimHitsH,
											       const edm::Handle<CrossingFrame<PSimHit>> &etlSimHitsH,
											       mtd::MTDGeomUtil &geomTools)
  : productGetter_(&productGetter), btlSimHitsH_(btlSimHitsH), etlSimHitsH_(etlSimHitsH), geomTools_(geomTools) {}

//
//---member functions
//

reco::SimToTPCollectionMtd MtdSimLayerClusterToTPAssociatorByTrackIdImpl::associateSimToTP(
      const edm::Handle<MtdSimLayerClusterCollection>& simClusH,
      const edm::Handle<TrackingParticleCollection>& trackingParticleH) const {

  
  SimToTPCollectionMtd outputCollection(productGetter_);

  // -- Fill maps of trackIdOffset per cell 
  std::unordered_map<uint64_t, std::vector<int> > map_simHitTrackId;
  std::unordered_map<uint64_t, std::vector<int> > map_simHitOffsetTrackId;
  fillSimHitTkIdMap(btlSimHitsH_,etlSimHitsH_, geomTools_, map_simHitTrackId, map_simHitOffsetTrackId);

  // -- get the collections
  const auto& simClusters  = *simClusH.product();
  const auto& trackingParticles = *trackingParticleH.product();

  // -- loop over sim clusters and get the trackId, eventId
  LogTrace("MtdSimLayerClusterToTPAssociator") << " Found " << simClusters.size() << " sim layer clusters in the event";
  
  size_t simClusIndex(0);
  for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
    auto simClus = *simClusIt;
    MtdSimLayerClusterRef simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex);
    unsigned int simClusTrackId = simClus.g4Tracks()[0].trackId();
    EncodedEventId simClusEventId = simClus.eventId();

    // -- Check the trackId offset of the sim hits and keep only clusters with "direct" hits (offset == 0)    
    /* to be added after implementing offset method in MtdSimLayerClusters
    /*
    std::vector<std::pair<uint64_t, float>> hitsAndFrac = simClus.hits_and_fractions();
    std::vector<uint64_t> simClusHitIds(hitsAndFrac.size());
    std::transform(hitsAndFrac.begin(), hitsAndFrac.end(), simClusHitIds.begin(), [](const std::pair<uint64_t, float>& pair) {
										    return pair.first;});

    //std::cout << " Sim layer cluster " << simClusIndex <<  " has " << simClusHitIds.size() << " sim Hits" << std::endl;
    //for (unsigned int i = 0; i < simClusHitIds.size(); i++){
    //  std::cout << "   Sim clus hit Id = " << simClusHitIds[i] <<std::endl;
    //  for (auto offset : map_simHitOffsetTrackId[simClusHitIds[i]]){
    //  	std::cout << "      -- TrackIdOffset = " << offset  <<std::endl;
    //  }
    //}
    
    bool hasOffset = false;
    for (unsigned int i = 0; i < simClusHitIds.size(); i++){
      auto it = std::find_if(map_simHitOffsetTrackId[simClusHitIds[i]].begin(), map_simHitOffsetTrackId[simClusHitIds[i]].end(), [](int offset) { return offset != 0;});
      if (it != map_simHitOffsetTrackId[simClusHitIds[i]].end()) hasOffset = true;
    }
    if (hasOffset) continue;
    */
    
    // -- loop over the tracking particles
    size_t tpIndex(0);
    for (auto tpIt = trackingParticles.begin(); tpIt != trackingParticles.end(); tpIt++) {
      auto tp = *tpIt;
      TrackingParticleRef tpRef = edm::Ref<TrackingParticleCollection>(trackingParticleH, tpIndex);
      unsigned int tpTrackId = tp.g4Tracks()[0].trackId();
      EncodedEventId tpEventId = tp.eventId();

      // -- fill the output collection
      if ( tpTrackId == simClusTrackId && tpEventId.rawId() == simClusEventId.rawId() ){
	LogDebug("MtdSimLayerClusterToTPAssociator") << "MtdSimLayerCluster: index = " << simClusIndex << "   simClusTrackId = " << simClusTrackId << " simClusEventId = " << simClusEventId.rawId() <<  " simClusEta = "<< simClus.eta();
	LogDebug("MtdSimLayerClusterToTPAssociator") << "  --> Found associated tracking particle:  index = " << tpIndex << "    tpTrackId = " << tpTrackId << " tpEventId = " << tpEventId.rawId() <<std::endl;
	outputCollection.insert(simClusterRef, tpRef);
      }

      tpIndex++;
      
    }// -- end loop over tracking particles

    simClusIndex++;

  }// -- end loop over sim clus
      
  return outputCollection;
  
}



reco::TPToSimCollectionMtd MtdSimLayerClusterToTPAssociatorByTrackIdImpl::associateTPToSim(
      const edm::Handle<MtdSimLayerClusterCollection>& simClusH,
      const edm::Handle<TrackingParticleCollection>& trackingParticleH) const {

  TPToSimCollectionMtd outputCollection(productGetter_);

  // -- Fill maps of trackIdOffset per cell
  std::unordered_map<uint64_t, std::vector<int> > map_simHitTrackId;
  std::unordered_map<uint64_t, std::vector<int> > map_simHitOffsetTrackId;
  fillSimHitTkIdMap(btlSimHitsH_,etlSimHitsH_, geomTools_, map_simHitTrackId, map_simHitOffsetTrackId);

  // -- get the collections
  const auto& simClusters  = *simClusH.product();
  const auto& trackingParticles = *trackingParticleH.product();

  // -- Loop over the tracking particles
  size_t tpIndex(0);
  for (auto tpIt = trackingParticles.begin(); tpIt != trackingParticles.end(); tpIt++) {
    auto tp = *tpIt;
    TrackingParticleRef tpRef = edm::Ref<TrackingParticleCollection>(trackingParticleH, tpIndex);
    unsigned int tpTrackId = tp.g4Tracks()[0].trackId();
    EncodedEventId tpEventId = tp.eventId();

    // -- loop over MtdSimLayerClusters
    size_t simClusIndex(0);
    for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
      auto simClus = *simClusIt;
      MtdSimLayerClusterRef simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex);
      unsigned int simClusTrackId = simClus.g4Tracks()[0].trackId();
      EncodedEventId simClusEventId = simClus.eventId();
      
      // -- Check the trackId offset of the sim hits and keep only clusters with "direct" hits (offset == 0)    
      /* to be added after implementing offset method in MtdSimLayerClusters
      /*
      std::vector<std::pair<uint64_t, float>> hitsAndFrac = simClus.hits_and_fractions();
      std::vector<uint64_t> simClusHitIds(hitsAndFrac.size());
      std::transform(hitsAndFrac.begin(), hitsAndFrac.end(), simClusHitIds.begin(), [](const std::pair<uint64_t, float>& pair) {
										    return pair.first;});

      bool hasOffset = false;
      for (unsigned int i = 0; i < simClusHitIds.size(); i++){
	auto it = std::find_if(map_simHitOffsetTrackId[simClusHitIds[i]].begin(), map_simHitOffsetTrackId[simClusHitIds[i]].end(), [](int offset) { return offset != 0;});
	if (it != map_simHitOffsetTrackId[simClusHitIds[i]].end()) hasOffset = true;
      }
      if (hasOffset) continue;
      */
      
      // -- Fill the output collection
      if ( tpTrackId == simClusTrackId && tpEventId.rawId() == simClusEventId.rawId() ){
        LogDebug("MtdSimLayerClusterToTPAssociator") << "Tracking particle:  index = " << tpIndex << "  tpTrackId = " << tpTrackId << "  tpEventId = " << tpEventId.rawId();
	LogDebug("MtdSimLayerClusterToTPAssociator") << "--> Found associated MtdSimLayerCluster: index = " << simClusIndex << "   simClusTrackId = " << simClusTrackId << " simClusEventId = " << simClusEventId.rawId() <<  " simClusEta = "<< simClus.eta();

	outputCollection.insert(tpRef, simClusterRef);
      }
      simClusIndex++;
    
    } // -- end loop over sim clusters

    tpIndex++;
    
  }
  
  return outputCollection;

}
