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
  std::unordered_map<uint64_t, std::vector<int> > map_simHitOffsetTrackId;
  fillSimHitTkIdOffsetMap(btlSimHitsH_,etlSimHitsH_, geomTools_, map_simHitOffsetTrackId);

  // -- get the collections
  const auto& simClusters  = *simClusH.product();
  const auto& trackingParticles = *trackingParticleH.product();

  // -- loop over sim clusters and get the trackId, eventId
  std::cout << " *** Found " << simClusters.size() << " sim layer clusters in the event" << std::endl;
  
  size_t simClusIndex(0);
  for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
    auto simClus = *simClusIt;
    MtdSimLayerClusterRef simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex);
    unsigned int simClusTrackId = simClus.g4Tracks()[0].trackId();
    EncodedEventId simClusEventId = simClus.eventId();

    //std::cout << "sim cluster " << simClusIndex << "   simClusTrackId = " << simClusTrackId << " simClusEventId = " << simClusEventId.rawId() <<std::endl;

    // -- Check the trackId offset of the sim hits and keep only clusters with "direct" hits (offset == 0)    
    std::vector<std::pair<uint64_t, float>> hitsAndFrac = simClus.hits_and_fractions();
    std::vector<uint64_t> simClusHitIds(hitsAndFrac.size());
    std::transform(hitsAndFrac.begin(), hitsAndFrac.end(), simClusHitIds.begin(), [](const std::pair<uint64_t, float>& pair) {
										    return pair.first;});

    //for (unsigned int i = 0; i < simClusHitIds.size(); i++){
    //  std::cout << "   Sim clus hit Id = " << simClusHitIds[i] <<std::endl;
    //  for (auto offset : map_simHitOffsetTrackId[simClusHitIds[i]]){
    //  	std::cout << "      -- TrackIdOffset = " << offset  <<std::endl;
    //  }
    // }
    
    bool hasOffset = false;
    for (unsigned int i = 0; i < simClusHitIds.size(); i++){
      auto it = std::find_if(map_simHitOffsetTrackId[simClusHitIds[i]].begin(), map_simHitOffsetTrackId[simClusHitIds[i]].end(), [](int offset) { return offset != 0;});
      if (it != map_simHitOffsetTrackId[simClusHitIds[i]].end()) hasOffset = true;
    }
    
    if (hasOffset) continue;
        
    // -- loop over the tracking particles
    size_t tpIndex(0);
    for (auto tpIt = trackingParticles.begin(); tpIt != trackingParticles.end(); tpIt++) {
      auto tp = *tpIt;
      TrackingParticleRef tpRef = edm::Ref<TrackingParticleCollection>(trackingParticleH, tpIndex);
      unsigned int tpTrackId = tp.g4Tracks()[0].trackId();
      EncodedEventId tpEventId = tp.eventId();

      // -- fill the output collection
      if ( tpTrackId == simClusTrackId && tpEventId.rawId() == simClusEventId.rawId() ){
	std::cout << " MTD SimLayerCluster " << simClusIndex << "   simClusTrackId = " << simClusTrackId << " simClusEventId = " << simClusEventId.rawId() <<  " simClusEta = "<< simClus.eta() <<  " has trackIdOffset = "<< hasOffset <<std::endl;
	std::cout << " Found associated tracking particle:  index = " << tpIndex << "    tpTrackId = " << tpTrackId << " tpEventId = " << tpEventId.rawId() <<std::endl;
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
  std::unordered_map<uint64_t, std::vector<int> > map_simHitOffsetTrackId;
  fillSimHitTkIdOffsetMap(btlSimHitsH_,etlSimHitsH_, geomTools_, map_simHitOffsetTrackId);

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

    //std::cout << "*** tracking particle " << tpIndex << "   tpTrackId = " << tpTrackId << " tpEventId = " << tpEventId.rawId() <<std::endl;

    // -- loop over MtdSimLayerClusters
    size_t simClusIndex(0);
    for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
      auto simClus = *simClusIt;
      MtdSimLayerClusterRef simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex);
      unsigned int simClusTrackId = simClus.g4Tracks()[0].trackId();
      EncodedEventId simClusEventId = simClus.eventId();
      
      // -- Check the trackId offset of the sim hits and keep only clusters with "direct" hits (offset == 0)    
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
        
      // -- Fill the output collection
      if ( tpTrackId == simClusTrackId && tpEventId.rawId() == simClusEventId.rawId() ){
	outputCollection.insert(tpRef, simClusterRef);
      }
      simClusIndex++;
    
    } // -- end loop over sim clusters

    tpIndex++;
    
  }
  
  return outputCollection;

}
