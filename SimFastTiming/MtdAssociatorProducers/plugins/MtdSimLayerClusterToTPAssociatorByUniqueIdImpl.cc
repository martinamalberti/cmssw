//
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "MtdSimLayerClusterToTPAssociatorByUniqueIdImpl.h"


using namespace reco;
using namespace std;

//using uniqueId = std::pair<unsigned int, EncodedEventId>;

/* Constructor */

MtdSimLayerClusterToTPAssociatorByUniqueIdImpl::MtdSimLayerClusterToTPAssociatorByUniqueIdImpl(edm::EDProductGetter const& productGetter)
  : productGetter_(&productGetter){}


//
//---member functions
//

reco::SimToTPCollectionMtd MtdSimLayerClusterToTPAssociatorByUniqueIdImpl::associateSimToTP(
      const edm::Handle<MtdSimLayerClusterCollection>& simClusH,
      const edm::Handle<TrackingParticleCollection>& trackingParticleH) const {

  SimToTPCollectionMtd outputCollection(productGetter_);
  
  // -- get the collections
  const auto& simClusters  = *simClusH.product();
  const auto& trackingParticles = *trackingParticleH.product();
  
  // -- loop over sim clusters and get the trackId, eventId
  size_t simClusIndex(0);
  for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
    auto simClus = *simClusIt;
    MtdSimLayerClusterRef simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex);
    unsigned int simClusTrackId = simClus.g4Tracks()[0].trackId();
    EncodedEventId simClusEventId = simClus.eventId();

    std::cout << "sim cluster " << simClusIndex << "   simClusTrackId = " << simClusTrackId << " simClusEventId = " << simClusEventId.rawId() <<std::endl;

    // check also the tkId offset of the sim hit and keep only clusters with "direct" hits (offset == 0)
    
    // -- loop over the tracking particles
    size_t tpIndex(0);
    for (auto tpIt = trackingParticles.begin(); tpIt != trackingParticles.end(); tpIt++) {
      auto tp = *tpIt;
      TrackingParticleRef tpRef = edm::Ref<TrackingParticleCollection>(trackingParticleH, tpIndex);
      unsigned int tpTrackId = tp.g4Tracks()[0].trackId();
      EncodedEventId tpEventId = tp.eventId();

      std::cout << "tracking particle " << tpIndex << "   tpTrackId = " << tpTrackId << " tpEventId = " << tpEventId.rawId() <<std::endl;

      // -- fill the output collection
      if ( tpTrackId == simClusTrackId && tpEventId.rawId() == simClusEventId.rawId() ){
	outputCollection.insert(simClusterRef, tpRef);
      }

      tpIndex++;
      
    }// -- end loop over tracking particles

    simClusIndex++;

  }// -- end loop over sim clus
      
  return outputCollection;
  
}



reco::TPToSimCollectionMtd MtdSimLayerClusterToTPAssociatorByUniqueIdImpl::associateTPToSim(
      const edm::Handle<MtdSimLayerClusterCollection>& simClusH,
      const edm::Handle<TrackingParticleCollection>& trackingParticleH) const {

  TPToSimCollectionMtd outputCollection(productGetter_);
  
  
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

    std::cout << "*** tracking particle " << tpIndex << "   tpTrackId = " << tpTrackId << " tpEventId = " << tpEventId.rawId() <<std::endl;

    // -- loop over MtdSimLayerClusters
    size_t simClusIndex(0);
    for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
      auto simClus = *simClusIt;
      MtdSimLayerClusterRef simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex);
      unsigned int simClusTrackId = simClus.g4Tracks()[0].trackId();
      EncodedEventId simClusEventId = simClus.eventId();

      // check also the tkId offset of the sim hit and keep only clusters with "direct" hits (offset == 0)
      
      // -- Fill th eouput collection
      if ( tpTrackId == simClusTrackId && tpEventId.rawId() == simClusEventId.rawId() ){
	outputCollection.insert(tpRef, simClusterRef);
      }
      simClusIndex++;
    
    } // -- end loop over sim clusters

    tpIndex++;
    
  }
  
  return outputCollection;

}
