//
//
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "SimFastTiming/MtdAssociatorProducers/interface/MtdRecoClusterToSimLayerClusterAssociatorImpl.h"

using namespace reco;
using namespace std;


/* Constructor */

MtdRecoClusterToSimLayerClusterAssociatorImpl(edm::EDProductGetter const& productGetter)
  : productGetter_(&productGetter){}


//
//---member functions
//

reco:RecoToSimCollection MtdRecoClusterToSimLayerClusterAssociatorImpl::associateRecoToSim(
  const edm::Handle<FTLlusterCollection>& rCCH, const edm::Handle<MtdSimLayerClusterCollection>& sCCH) const {

  RecoToSimCollection outputCollection(productGetter_);

  // do stuff to fill the output collection
  

  return outputCollection;
}



reco:SimToRecoCollection MtdRecoClusterToSimLayerClusterAssociatorImpl::associateSimToReco(
  const edm::Handle<FTLlusterCollection>& rCCH, const edm::Handle<MtdSimLayerClusterCollection>& sCCH) const {

  SimToRecoCollection outputCollection(productGetter_);

  // do stuff to fill the output collection
  


  return outputCollection;
}
