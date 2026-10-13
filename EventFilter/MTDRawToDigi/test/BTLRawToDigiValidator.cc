// BTLRawToDigiValidator.cc
//
// EDAnalyzer that compares, event by event and digi by digi:
//   - refDigis  : BTLDigiContentCollection (edm::SortedCollection<btldigi::BTLDigi>)
//                 produced by the digitizer (mix, MTDBarrel)
//   - testDigis : btldigi::BTLDigiHostCollection (PortableHostCollection<BTLDigiSoA>)
//                 produced by BTLRawToDigi (unpacker output)
//
// Assumes both collections are in the same order (same index = same digi).
//

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/FTLDigi/interface/MTDDigiCollections.h"

// SoA host collection from the unpacker
#include "DataFormats/FTLDigiSoA/interface/BTLDigiSoA.h"
#include "DataFormats/FTLDigiSoA/interface/BTLDigiHostCollection.h"

class BTLRawToDigiValidator : public edm::one::EDAnalyzer<> {
public:
  explicit BTLRawToDigiValidator(const edm::ParameterSet&);
  ~BTLRawToDigiValidator() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ---------- member data ----------
  edm::EDGetTokenT<BTLDigiContentCollection>       refDigisToken_;   // edm::SortedCollection<btldigi::BTLDigi>
  edm::EDGetTokenT<btldigi::BTLDigiHostCollection> testDigisToken_;  // SoA host collection
};

// ---------------------------------------------------------------------------
// Constructor / destructor
// ---------------------------------------------------------------------------

BTLRawToDigiValidator::BTLRawToDigiValidator(const edm::ParameterSet& iConfig) {
  refDigisToken_  = consumes<BTLDigiContentCollection>(iConfig.getParameter<edm::InputTag>("refDigis"));
  testDigisToken_ = consumes<btldigi::BTLDigiHostCollection>(iConfig.getParameter<edm::InputTag>("testDigis"));
}

BTLRawToDigiValidator::~BTLRawToDigiValidator() {}

// ---------------------------------------------------------------------------
// beginJob / endJob
// ---------------------------------------------------------------------------

void BTLRawToDigiValidator::beginJob() {}

void BTLRawToDigiValidator::endJob() {}

// ---------------------------------------------------------------------------
// analyze
// ---------------------------------------------------------------------------

void BTLRawToDigiValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<BTLDigiContentCollection>       h_ref;
  edm::Handle<btldigi::BTLDigiHostCollection> h_test;

  iEvent.getByToken(refDigisToken_,  h_ref);
  iEvent.getByToken(testDigisToken_, h_test);

  auto const testView = h_test->view();  // BTLDigiSoA::ConstView

  int32_t const nRef  = static_cast<int32_t>(h_ref->size());
  int32_t const nTest = testView.metadata().size();

  edm::LogInfo("BTLRawToDigiValidator")
    << "Event " << iEvent.id() << "  nDigis(ref)=" << nRef << "  nDigis(test)=" << nTest;

  if (nRef != nTest) {
    edm::LogError("BTLRawToDigiValidator")
      << "Event " << iEvent.id() << "  SIZE MISMATCH: ref=" << nRef << " test=" << nTest;
  }

  int32_t const nCompare = std::min(nRef, nTest);
  int nOk = 0, nWrong = 0;

  for (int32_t i = 0; i < nCompare; ++i) {
    btldigi::BTLDigi const& r = (*h_ref)[i];
    auto const t = testView[i];
    
    bool ok = true;  // reset for each digi
    
    // Lambda: compares one field by value and logs an error on mismatch.
    auto check = [&](const char* fieldName, auto refVal, auto testVal) {
		   if (refVal != testVal) {
		     edm::LogError("BTLRawToDigiValidator")
		     << "Event " << iEvent.id() << "   digi[" << i << "]"
		     << "  mismatch in " << fieldName
		     << ": ref=" << +refVal << " test=" << +testVal;
		     ok = false;
		   }
		 };
    

    check("rawId",         r.krawId(),         t.rawId());
    check("BC0count",      r.kBC0count(),      t.BC0count());
    check("status",        r.kstatus(),        t.status());
    check("BCcount",       r.kBCcount(),       t.BCcount());
    // plus side
    check("chIDPlus",      r.kchIDPlus(),      t.chIDPlus());
    check("T1coarsePlus",  r.kT1coarsePlus(),  t.T1coarsePlus());
    check("T2coarsePlus",  r.kT2coarsePlus(),  t.T2coarsePlus());
    check("EOIcoarsePlus", r.kEOIcoarsePlus(), t.EOIcoarsePlus());
    check("ChargePlus",    r.kChargePlus(),    t.ChargePlus());
    check("T1finePlus",    r.kT1finePlus(),    t.T1finePlus());
    check("T2finePlus",    r.kT2finePlus(),    t.T2finePlus());
    check("IdleTimePlus",  r.kIdleTimePlus(),  t.IdleTimePlus());
    check("PrevTrigFPlus", r.kPrevTrigFPlus(), t.PrevTrigFPlus());
    check("TACIDPlus",     r.kTACIDPlus(),     t.TACIDPlus());
    // minus side
    check("chIDMinus",      r.kchIDMinus(),      t.chIDMinus());
    check("T1coarseMinus",  r.kT1coarseMinus(),  t.T1coarseMinus());
    check("T2coarseMinus",  r.kT2coarseMinus(),  t.T2coarseMinus());
    check("EOIcoarseMinus", r.kEOIcoarseMinus(), t.EOIcoarseMinus());
    check("ChargeMinus",    r.kChargeMinus(),    t.ChargeMinus());
    check("T1fineMinus",    r.kT1fineMinus(),    t.T1fineMinus());
    check("T2fineMinus",    r.kT2fineMinus(),    t.T2fineMinus());
    check("IdleTimeMinus",  r.kIdleTimeMinus(),  t.IdleTimeMinus());
    check("PrevTrigFMinus", r.kPrevTrigFMinus(), t.PrevTrigFMinus());
    check("TACIDMinus",     r.kTACIDMinus(),     t.TACIDMinus());
    
    if (ok) ++nOk;
    else ++nWrong;
  }
  
  if (nRef == nTest && nWrong == 0)
    edm::LogInfo("BTLRawToDigiValidator")
      << "Event " << iEvent.id() << "  OK - all " << nOk << " digis match.";
  else
    edm::LogError("BTLRawToDigiValidator")
      << "Event " << iEvent.id() << " FAILURES:"
      << "  matched=" << nOk << "  wrong=" << nWrong
      << "  size_diff=" << (nRef - nTest);
}

// ---------------------------------------------------------------------------
// fillDescriptions
// ---------------------------------------------------------------------------

void BTLRawToDigiValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("refDigis",  edm::InputTag{"mix", "MTDBarrel"});
  desc.add<edm::InputTag>("testDigis", edm::InputTag{"btlRawToDigiGPU"});
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(BTLRawToDigiValidator);
