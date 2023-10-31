#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/ValidHandle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/GeantUnits.h"
#include "DataFormats/Math/interface/angle_units.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"
#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/DetLayers/interface/MTDTrayBarrelLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetTray.h"
#include "RecoMTD/DetLayers/interface/MTDSectorForwardDoubleLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetSector.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"

#include "DataFormats/Common/interface/OneToMany.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "HepMC/GenRanges.h"
#include "CLHEP/Units/PhysicalConstants.h"


#include "MTDHit.h"
#include "MTDTrack.h"

#include "TTree.h"



class MtdTracksDumper : public edm::one::EDAnalyzer<edm::one::WatchRuns,
						    edm::one::SharedResources> {
public:
  explicit MtdTracksDumper(const edm::ParameterSet&);
  ~MtdTracksDumper() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override{};
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override{};
  void endJob() override ;
    
  void clearTreeVariables();

  const bool mvaGenSel(const HepMC::GenParticle&, const float&);
  const bool mvaTPSel(const TrackingParticle&);
  const bool mvaRecSel(const reco::TrackBase&, const reco::Vertex&, const double&, const double&);
  const bool mvaGenRecMatch(const HepMC::GenParticle&, const double&, const reco::TrackBase&, const bool&);
  const edm::Ref<std::vector<TrackingParticle>>* getMatchedTP(const reco::TrackBaseRef&);

  // ------------ member data ------------
  const float trackMinPt_;
  const float trackMaxBtlEta_;
  const float trackMinEtlEta_;
  const float trackMaxEtlEta_;

  static constexpr double etacutGEN_ = 4.;               // |eta| < 4;
  static constexpr double etacutREC_ = 3.;               // |eta| < 3;
  static constexpr double pTcut_ = 0.7;                  // PT > 0.7 GeV
  static constexpr double deltaZcut_ = 0.1;              // dz separation 1 mm
  static constexpr double deltaPTcut_ = 0.05;            // dPT < 5%
  static constexpr double deltaDRcut_ = 0.03;            // DeltaR separation
  static constexpr double rBTL_ = 110.0;
  static constexpr double zETL_ = 290.0;
  
  const reco::RecoToSimCollection* r2s_;
  const reco::SimToRecoCollection* s2r_;

  edm::EDGetTokenT<reco::TrackCollection> GenRecTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> RecTrackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> RecVertexToken_;

  edm::EDGetTokenT<edm::HepMCProduct> HepMCProductToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoAssociationToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;

  edm::EDGetTokenT<edm::ValueMap<int>> trackAssocToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> pathLengthToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmatmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0PidToken_;
  //edm::EDGetTokenT<edm::ValueMap<float>> t0SafePidToken_;
  //edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> trackMVAQualToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchTimeChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> etlMatchChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> etlMatchTimeChi2Token_;
    
  edm::ESGetToken<MTDGeometry, MTDDigiGeometryRecord> mtdgeoToken_;
  edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdtopoToken_;
  edm::ESGetToken<MTDDetLayerGeometry, MTDRecoGeometryRecord> mtdlayerToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magfieldToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> builderToken_;
  edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particleTableToken_;

  // --  Output tree
  TTree* tree_; 
  MTDTrack *mtdTrackInfo;  

};

// ------------ constructor and destructor --------------
MtdTracksDumper::MtdTracksDumper(const edm::ParameterSet& iConfig):
  trackMinPt_(iConfig.getParameter<double>("trackMinimumPt")),
  trackMaxBtlEta_(iConfig.getParameter<double>("trackMaximumBtlEta")),
  trackMinEtlEta_(iConfig.getParameter<double>("trackMinimumEtlEta")),
  trackMaxEtlEta_(iConfig.getParameter<double>("trackMaximumEtlEta"))
{
  GenRecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagG"));
  RecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagT"));
  RecVertexToken_ = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("inputTagV"));
  HepMCProductToken_ = consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("inputTagH"));
  trackingParticleCollectionToken_ = consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  simToRecoAssociationToken_ = consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  recoToSimAssociationToken_ = consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  trackAssocToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackAssocSrc"));
  pathLengthToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pathLengthSrc"));
  tmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtd"));
  SigmatmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatmtd"));
  t0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"));
  Sigmat0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"));
  t0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0PID"));
  Sigmat0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0PID"));
  //t0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0SafePID"));
  // Sigmat0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0SafePID"));
  trackMVAQualToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMVAQual"));
  btlMatchChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchChi2Tag"));
  btlMatchTimeChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchTimeChi2Tag"));
  etlMatchChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchChi2Tag"));
  etlMatchTimeChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchTimeChi2Tag"));
  mtdgeoToken_ = esConsumes<MTDGeometry, MTDDigiGeometryRecord>();
  mtdtopoToken_ = esConsumes<MTDTopology, MTDTopologyRcd>();
  mtdlayerToken_ = esConsumes<MTDDetLayerGeometry, MTDRecoGeometryRecord>();
  magfieldToken_ = esConsumes<MagneticField, IdealMagneticFieldRecord>();
  builderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));
  particleTableToken_ = esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>();

}

// ------------ method called once each job just before starting event loop  ------------
void MtdTracksDumper::beginJob() {

  mtdTrackInfo = new MTDTrack;
  
  edm::Service<TFileService> fs_;
  tree_ = fs_->make<TTree>("tracks", "tracks-MTD tree");

  tree_ -> Branch( "pv_x",       &mtdTrackInfo->pv_x);
  tree_ -> Branch( "pv_y",       &mtdTrackInfo->pv_y);
  tree_ -> Branch( "pv_z",       &mtdTrackInfo->pv_z);
  tree_ -> Branch( "pv_t",       &mtdTrackInfo->pv_t);
  tree_ -> Branch( "pv_isFake",  &mtdTrackInfo->pv_isFake);

  tree_ -> Branch( "genpv_x",    &mtdTrackInfo->genpv_x);
  tree_ -> Branch( "genpv_y",    &mtdTrackInfo->genpv_y);
  tree_ -> Branch( "genpv_z",    &mtdTrackInfo->genpv_z);
  tree_ -> Branch( "genpv_t",    &mtdTrackInfo->genpv_t);

  tree_ -> Branch( "p",          &mtdTrackInfo->p);
  tree_ -> Branch( "pt",         &mtdTrackInfo->pt);
  tree_ -> Branch( "eta",        &mtdTrackInfo->eta);
  tree_ -> Branch( "phi",        &mtdTrackInfo->phi);
  tree_ -> Branch( "x",          &mtdTrackInfo->x);
  tree_ -> Branch( "y",          &mtdTrackInfo->y);
  tree_ -> Branch( "z",          &mtdTrackInfo->z);
  tree_ -> Branch( "dz",         &mtdTrackInfo->dz);
  tree_ -> Branch( "dxy",        &mtdTrackInfo->dxy);
  tree_ -> Branch( "dzErr",      &mtdTrackInfo->dzErr);
  tree_ -> Branch( "dxyErr",     &mtdTrackInfo->dxyErr);
  tree_ -> Branch( "chi2",       &mtdTrackInfo->chi2);
  tree_ -> Branch( "ndof",       &mtdTrackInfo->ndof);
  tree_ -> Branch( "nTrackerHits",       &mtdTrackInfo->nTrackerHits);
  tree_ -> Branch( "nPixelBarrelHits",   &mtdTrackInfo->nPixelBarrelHits);
  tree_ -> Branch( "nPixelEndcapHits",   &mtdTrackInfo->nPixelEndcapHits);  
  tree_ -> Branch( "nStripHits",         &mtdTrackInfo->nStripHits);
  tree_ -> Branch( "nTimingBTLHits",     &mtdTrackInfo->nTimingBTLHits);
  tree_ -> Branch( "nTimingETLHits",     &mtdTrackInfo->nTimingETLHits);

  tree_ -> Branch( "btlMatchChi2",       &mtdTrackInfo->btlMatchChi2);
  tree_ -> Branch( "btlMatchTimeChi2",   &mtdTrackInfo->btlMatchTimeChi2);
  tree_ -> Branch( "etlMatchChi2",       &mtdTrackInfo->etlMatchChi2);
  tree_ -> Branch( "etlMatchTimeChi2",   &mtdTrackInfo->etlMatchTimeChi2);
  
  tree_ -> Branch( "tmtd",               &mtdTrackInfo->tmtd);
  tree_ -> Branch( "t0",                 &mtdTrackInfo->t0);
  tree_ -> Branch( "sigmat0",            &mtdTrackInfo->sigmat0);
  tree_ -> Branch( "t0Pid",              &mtdTrackInfo->t0Pid);
  tree_ -> Branch( "sigmat0Pid",         &mtdTrackInfo->sigmat0Pid);
  tree_ -> Branch( "pathlength",         &mtdTrackInfo->pathlength);

  tree_ -> Branch( "mtdTrackQualityMVA", &mtdTrackInfo->mtdTrackQualityMVA);

  tree_ -> Branch( "isMatchedToTP",      &mtdTrackInfo->isMatchedToTP);
  tree_ -> Branch( "isMatchedToSimCluster", &mtdTrackInfo->isMatchedToSimCluster);
  tree_ -> Branch( "isMatchedToGen",       &mtdTrackInfo->isMatchedToGen);
  

}


MtdTracksDumper::~MtdTracksDumper() {}


// ------------ method called for each event  ------------
void MtdTracksDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // -- initialize output tree
  clearTreeVariables();

  using namespace edm;
  using namespace geant_units::operators;
  using namespace std;

  // -- General tracks 
  auto GenRecTrackHandle = makeValid(iEvent.getHandle(GenRecTrackToken_));

  // -- reco to tracking particles associations
  auto simToRecoH = makeValid(iEvent.getHandle(simToRecoAssociationToken_));
  s2r_ = simToRecoH.product();

  auto recoToSimH = makeValid(iEvent.getHandle(recoToSimAssociationToken_));
  r2s_ = recoToSimH.product();

  // -- MTD related infos
  const auto& tMtd = iEvent.get(tmtdToken_);
  const auto& SigmatMtd = iEvent.get(SigmatmtdToken_);
  const auto& t0Src = iEvent.get(t0SrcToken_);
  const auto& Sigmat0Src = iEvent.get(Sigmat0SrcToken_);
  const auto& t0PidMap = iEvent.get(t0PidToken_);
  const auto& Sigmat0PidMap = iEvent.get(Sigmat0PidToken_);
  //const auto& t0Safe = iEvent.get(t0SafePidToken_);
  //const auto& Sigmat0Safe = iEvent.get(Sigmat0SafePidToken_);
  const auto& mtdQualMVA = iEvent.get(trackMVAQualToken_);
  const auto& trackAssoc = iEvent.get(trackAssocToken_);
  const auto& pathLength = iEvent.get(pathLengthToken_);
  const auto& btlMatchChi2Map = iEvent.get(btlMatchChi2Token_);
  const auto& btlMatchTimeChi2Map = iEvent.get(btlMatchTimeChi2Token_);
  const auto& etlMatchChi2Map = iEvent.get(etlMatchChi2Token_);
  const auto& etlMatchTimeChi2Map = iEvent.get(etlMatchTimeChi2Token_);
  
  // -- Reco vertices
  auto RecVertexHandle = makeValid(iEvent.getHandle(RecVertexToken_));
  auto recoVtx = *RecVertexHandle.product();
   
  // -- generator level information (HepMC format)
  auto GenEventHandle = makeValid(iEvent.getHandle(HepMCProductToken_));
  const HepMC::GenEvent* mc = GenEventHandle->GetEvent();
  double xsim = convertMmToCm((*(mc->vertices_begin()))->position().x());
  double ysim = convertMmToCm((*(mc->vertices_begin()))->position().y());
  double zsim = convertMmToCm((*(mc->vertices_begin()))->position().z());
  double tsim = (*(mc->vertices_begin()))->position().t() * CLHEP::mm / CLHEP::c_light;

  auto pdt = iSetup.getHandle(particleTableToken_);
  const HepPDT::ParticleDataTable* pdTable = pdt.product();


  // -- Fill tree with reco PV info - 4D reco vertices
  if (recoVtx.size() > 0){
    mtdTrackInfo->pv_x = recoVtx[0].x();
    mtdTrackInfo->pv_y = recoVtx[0].y();
    mtdTrackInfo->pv_z = recoVtx[0].z();
    mtdTrackInfo->pv_t = recoVtx[0].t();
    mtdTrackInfo->pv_isFake = recoVtx[0].isFake();
  }

  // -- Fill tree with gen vtx info
  mtdTrackInfo->genpv_x = xsim;
  mtdTrackInfo->genpv_y = ysim;
  mtdTrackInfo->genpv_z = zsim;
  mtdTrackInfo->genpv_t = tsim;
  
  unsigned int index = 0;  
  // --- Loop over all RECO tracks ---
  for (const auto& trackGen : *GenRecTrackHandle) {
    const reco::TrackRef trackref(iEvent.getHandle(GenRecTrackToken_), index);
    index++;

    if (trackAssoc[trackref] == -1) {
      LogInfo("mtdTracks") << "Extended track not associated";
      continue;
    }

    const reco::TrackRef mtdTrackref = reco::TrackRef(iEvent.getHandle(RecTrackToken_), trackAssoc[trackref]);
    const reco::Track& track = *mtdTrackref;

    if (track.pt() >= trackMinPt_ && std::abs(track.eta()) <= trackMaxEtlEta_) {

      // - Fill tree with track infos
      mtdTrackInfo->x.push_back(track.vx());
      mtdTrackInfo->y.push_back(track.vy());
      mtdTrackInfo->z.push_back(track.vz());
      mtdTrackInfo->dz.push_back(track.dz(recoVtx[0].position()));
      mtdTrackInfo->dxy.push_back(track.dxy(recoVtx[0].position()));
      mtdTrackInfo->dzErr.push_back(track.dzError());//? correct if PV is the first vtx in the 4D vertex collection?
      mtdTrackInfo->dxyErr.push_back(track.dxyError()); //? correct if PV is the first vtx in the 4D vertex collection?
      mtdTrackInfo->p.push_back(track.p());
      mtdTrackInfo->pt.push_back(track.pt());
      mtdTrackInfo->eta.push_back(track.eta());
      mtdTrackInfo->phi.push_back(track.phi());
      mtdTrackInfo->chi2.push_back(track.chi2());//?
      mtdTrackInfo->ndof.push_back(track.ndof());//?
      mtdTrackInfo->nTrackerHits.push_back(track.hitPattern().numberOfValidTrackerHits());
      mtdTrackInfo->nPixelBarrelHits.push_back(track.hitPattern().numberOfValidPixelBarrelHits());
      mtdTrackInfo->nPixelEndcapHits.push_back(track.hitPattern().numberOfValidPixelEndcapHits());
      mtdTrackInfo->nStripHits.push_back(track.hitPattern().numberOfValidStripHits());
      mtdTrackInfo->nTimingBTLHits.push_back(track.hitPattern().numberOfValidTimingBTLHits());
      mtdTrackInfo->nTimingETLHits.push_back(track.hitPattern().numberOfValidTimingETLHits());
      mtdTrackInfo->btlMatchChi2.push_back(btlMatchChi2Map[trackref]);
      mtdTrackInfo->btlMatchTimeChi2.push_back(btlMatchTimeChi2Map[trackref]);
      mtdTrackInfo->etlMatchChi2.push_back(etlMatchChi2Map[trackref]);
      mtdTrackInfo->etlMatchTimeChi2.push_back(etlMatchTimeChi2Map[trackref]);
      mtdTrackInfo->tmtd.push_back(tMtd[trackref]);
      mtdTrackInfo->t0.push_back(t0Src[trackref]);
      mtdTrackInfo->sigmat0.push_back(Sigmat0Src[trackref]);
      mtdTrackInfo->t0Pid.push_back(t0PidMap[trackref]);
      mtdTrackInfo->sigmat0Pid.push_back(Sigmat0PidMap[trackref]);
      mtdTrackInfo->pathlength.push_back(pathLength[trackref]);
      mtdTrackInfo->mtdTrackQualityMVA.push_back(mtdQualMVA[trackref]);

      // - TrackingParticle based matching
      const reco::TrackBaseRef trkrefb(trackref);
      auto tp_info = getMatchedTP(trkrefb);

      if (tp_info != nullptr && mvaTPSel(**tp_info) ) mtdTrackInfo->isMatchedToTP.push_back(1);
      else mtdTrackInfo->isMatchedToTP.push_back(0);

      // -- Matching with gen particles
      int genMatching = 0;
      for (const auto& genP : mc->particle_range()) {
          // select status 1 genParticles and match them to the reconstructed track
	float charge = pdTable->particle(HepPDT::ParticleID(genP->pdg_id())) != nullptr ? pdTable->particle(HepPDT::ParticleID(genP->pdg_id()))->charge() : 0.f;
	if (mvaGenSel(*genP, charge) &&  mvaGenRecMatch(*genP, zsim, trackGen, recoVtx[0].isFake())) genMatching = 1;
      }
      mtdTrackInfo->isMatchedToGen.push_back(genMatching);
 
      // -- matching with SimLayerCluster.. to DO!
      

    }
  } // -- end RECO tracks loop

  std::cout << "------------------------------\n";

  tree_->Fill();

}



void MtdTracksDumper::clearTreeVariables() {

  mtdTrackInfo->pv_x =  -999;
  mtdTrackInfo->pv_y =  -999;
  mtdTrackInfo->pv_z =  -999;
  mtdTrackInfo->pv_t =  -999;
  mtdTrackInfo->pv_isFake =  -999;

  mtdTrackInfo->genpv_x =  -999;
  mtdTrackInfo->genpv_y =  -999;
  mtdTrackInfo->genpv_z =  -999;
  mtdTrackInfo->genpv_t =  -999;
  
  mtdTrackInfo->x.clear();
  mtdTrackInfo->y.clear();
  mtdTrackInfo->z.clear();
  mtdTrackInfo->dz.clear();
  mtdTrackInfo->dxy.clear();
  mtdTrackInfo->dzErr.clear();
  mtdTrackInfo->dxyErr.clear();
  mtdTrackInfo->p.clear();
  mtdTrackInfo->pt.clear();
  mtdTrackInfo->eta.clear();
  mtdTrackInfo->phi.clear();
  mtdTrackInfo->chi2.clear();
  mtdTrackInfo->ndof.clear();
  mtdTrackInfo->nTrackerHits.clear();
  mtdTrackInfo->nPixelBarrelHits.clear();
  mtdTrackInfo->nPixelEndcapHits.clear();
  mtdTrackInfo->nStripHits.clear();
  mtdTrackInfo->nTimingBTLHits.clear();
  mtdTrackInfo->nTimingETLHits.clear();

  mtdTrackInfo->btlMatchChi2.clear();
  mtdTrackInfo->btlMatchTimeChi2.clear();
  mtdTrackInfo->etlMatchChi2.clear();
  mtdTrackInfo->etlMatchTimeChi2.clear();
  
  mtdTrackInfo->pathlength.clear();
  mtdTrackInfo->tmtd.clear();
  mtdTrackInfo->t0.clear();
  mtdTrackInfo->sigmat0.clear();
  mtdTrackInfo->t0Pid.clear();
  mtdTrackInfo->sigmat0Pid.clear();
  mtdTrackInfo->mtdTrackQualityMVA.clear();

  mtdTrackInfo->isMatchedToTP.clear();
  mtdTrackInfo->isMatchedToSimCluster.clear();
  mtdTrackInfo->isMatchedToGen.clear();

}
  



void MtdTracksDumper::endJob() {}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void MtdTracksDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/Tracks");
  desc.add<edm::InputTag>("inputTagG", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("inputTagT", edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>("inputTagV", edm::InputTag("offlinePrimaryVertices4D"));
  desc.add<edm::InputTag>("inputTagH", edm::InputTag("generatorSmeared"));
  desc.add<edm::InputTag>("TPtoRecoTrackAssoc", edm::InputTag("trackingParticleRecoTrackAsssociation"));
  desc.add<edm::InputTag>("tmtd", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmatmtd", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("t0Src", edm::InputTag("trackExtenderWithMTD:generalTrackt0"));
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("trackExtenderWithMTD:generalTracksigmat0"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"))
      ->setComment("Association between General and MTD Extended tracks");
  desc.add<edm::InputTag>("pathLengthSrc", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("t0SafePID", edm::InputTag("tofPID:t0safe"));
  desc.add<edm::InputTag>("sigmat0SafePID", edm::InputTag("tofPID:sigmat0safe"));
  desc.add<edm::InputTag>("sigmat0PID", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("t0PID", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("trackMVAQual", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<double>("trackMinimumPt", 0.7);  // [GeV]
  desc.add<double>("trackMaximumBtlEta", 1.5);
  desc.add<double>("trackMinimumEtlEta", 1.6);
  desc.add<double>("trackMaximumEtlEta", 3.);
}




const bool MtdTracksDumper::mvaGenSel(const HepMC::GenParticle& gp, const float& charge) {
  bool match = false;
  if (gp.status() != 1) {
    return match;
  }
  match = charge != 0.f && gp.momentum().perp() > pTcut_ && std::abs(gp.momentum().eta()) < etacutGEN_;
  return match;
}

const bool MtdTracksDumper::mvaTPSel(const TrackingParticle& tp) {
  bool match = false;
  if (tp.status() != 1) {
    return match;
  }
  auto x_pv = tp.parentVertex()->position().x();
  auto y_pv = tp.parentVertex()->position().y();
  auto z_pv = tp.parentVertex()->position().z();

  auto r_pv = std::sqrt(x_pv * x_pv + y_pv * y_pv);

  match = tp.charge() != 0 && tp.pt() > pTcut_ && std::abs(tp.eta()) < etacutGEN_ && r_pv < rBTL_ && z_pv < zETL_;
  return match;
}

const bool MtdTracksDumper::mvaRecSel(const reco::TrackBase& trk,
                                          const reco::Vertex& vtx,
                                          const double& t0,
                                          const double& st0) {
  bool match = false;
  match = trk.pt() > pTcut_ && std::abs(trk.eta()) < etacutREC_ &&
          (std::abs(trk.vz() - vtx.z()) <= deltaZcut_ || vtx.isFake());
  if (st0 > 0.) {
    match = match && std::abs(t0 - vtx.t()) < 3. * st0;
  }
  return match;
}

const bool MtdTracksDumper::mvaGenRecMatch(const HepMC::GenParticle& genP,
                                               const double& zsim,
                                               const reco::TrackBase& trk,
                                               const bool& vtxFake) {
  bool match = false;
  double dR = reco::deltaR(genP.momentum(), trk.momentum());
  double genPT = genP.momentum().perp();
  match = std::abs(genPT - trk.pt()) < trk.pt() * deltaPTcut_ && dR < deltaDRcut_ &&
          (std::abs(trk.vz() - zsim) < deltaZcut_ || vtxFake);
  return match;
}

const edm::Ref<std::vector<TrackingParticle>>* MtdTracksDumper::getMatchedTP(const reco::TrackBaseRef& recoTrack) {
  auto found = r2s_->find(recoTrack);

  // reco track not matched to any TP
  if (found == r2s_->end())
    return nullptr;

  //matched TP equal to any TP associated to in time events
  for (const auto& tp : found->val) {
    if (tp.first->eventId().bunchCrossing() == 0)
      return &tp.first;
  }

  // reco track not matched to any TP from vertex
  return nullptr;
}


DEFINE_FWK_MODULE(MtdTracksDumper);
