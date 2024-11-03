#include <string>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

class MtdTracksHarvester : public DQMEDHarvester {
public:
  explicit MtdTracksHarvester(const edm::ParameterSet& iConfig);
  ~MtdTracksHarvester() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:
  void dqmEndJob(DQMStore::IBooker&, DQMStore::IGetter&) override;

private:
  void computeEfficiency1D(MonitorElement* num, MonitorElement* den, MonitorElement* result);
  void normalize(MonitorElement* h, double scale);

  const std::string folder_;

  // --- Histograms
  MonitorElement* meBtlEtaEff_;
  MonitorElement* meBtlPhiEff_;
  MonitorElement* meBtlPtEff_;
  MonitorElement* meEtlEtaEff_;
  MonitorElement* meEtlPhiEff_;
  MonitorElement* meEtlPtEff_;
  MonitorElement* meEtlEtaEff2_;
  MonitorElement* meEtlPhiEff2_;
  MonitorElement* meEtlPtEff2_;
  MonitorElement* meEtlEtaEffLowPt_[2];
  MonitorElement* meEtlEtaEff2LowPt_[2];
  MonitorElement* meTPPtSelEff_;
  MonitorElement* meTPEtaSelEff_;
  MonitorElement* meTPPtMatchEff_;
  MonitorElement* meTPEtaMatchEff_;
  MonitorElement* meTPPtMatchEtl2Eff_;
  MonitorElement* meTPEtaMatchEtl2Eff_;

  // - BTL track-mtd matching efficiencies
  MonitorElement* meBTLTPmtdDirectEtaSelEff_; 
  MonitorElement* meBTLTPmtdDirectPtSelEff_;
  MonitorElement* meBTLTPmtdOtherEtaSelEff_; 
  MonitorElement* meBTLTPmtdOtherPtSelEff_;
  MonitorElement* meBTLTPnomtdEtaSelEff_; 
  MonitorElement* meBTLTPnomtdPtSelEff_;

  MonitorElement* meBTLTPmtdDirectCorrectAssocEtaMatchEff_; 
  MonitorElement* meBTLTPmtdDirectCorrectAssocPtMatchEff_;
  MonitorElement* meBTLTPmtdDirectWrongAssocEtaMatchEff_; 
  MonitorElement* meBTLTPmtdDirectWrongAssocPtMatchEff_;
  MonitorElement* meBTLTPmtdDirectNoAssocEtaMatchEff_; 
  MonitorElement* meBTLTPmtdDirectNoAssocPtMatchEff_;
  
  MonitorElement* meBTLTPmtdOtherCorrectAssocEtaMatchEff_; 
  MonitorElement* meBTLTPmtdOtherCorrectAssocPtMatchEff_;
  MonitorElement* meBTLTPmtdOtherWrongAssocEtaMatchEff_; 
  MonitorElement* meBTLTPmtdOtherWrongAssocPtMatchEff_;
  MonitorElement* meBTLTPmtdOtherNoAssocEtaMatchEff_; 
  MonitorElement* meBTLTPmtdOtherNoAssocPtMatchEff_;
  
  MonitorElement* meBTLTPnomtdEtaMatchEff_; 
  MonitorElement* meBTLTPnomtdPtMatchEff_;

  // - ETL track-mtd matching efficiencies
  MonitorElement* meETLTPmtd1EtaSelEff_;
  MonitorElement* meETLTPmtd1PtSelEff_;
  MonitorElement* meETLTPmtd2EtaSelEff_;
  MonitorElement* meETLTPmtd2PtSelEff_;
  MonitorElement* meETLTPnomtdEtaSelEff_;
  MonitorElement* meETLTPnomtdPtSelEff_;

  MonitorElement* meETLTPmtd1CorrectAssocEtaMatchEff_; 
  MonitorElement* meETLTPmtd1CorrectAssocPtMatchEff_;
  MonitorElement* meETLTPmtd1WrongAssocEtaMatchEff_; 
  MonitorElement* meETLTPmtd1WrongAssocPtMatchEff_;
  MonitorElement* meETLTPmtd1NoAssocEtaMatchEff_; 
  MonitorElement* meETLTPmtd1NoAssocPtMatchEff_;
  
  MonitorElement* meETLTPmtd2CorrectAssocEtaMatchEff_; 
  MonitorElement* meETLTPmtd2CorrectAssocPtMatchEff_;
  MonitorElement* meETLTPmtd2WrongAssocEtaMatchEff_; 
  MonitorElement* meETLTPmtd2WrongAssocPtMatchEff_;
  MonitorElement* meETLTPmtd2NoAssocEtaMatchEff_; 
  MonitorElement* meETLTPmtd2NoAssocPtMatchEff_;
  
  MonitorElement* meETLTPnomtdEtaMatchEff_; 
  MonitorElement* meETLTPnomtdPtMatchEff_;

  //
  MonitorElement* meNoTimeFraction_;
  MonitorElement* meExtraPtEff_;
  MonitorElement* meExtraPtEtl2Eff_;
  MonitorElement* meExtraEtaEff_;
  MonitorElement* meExtraEtaEtl2Eff_;
  MonitorElement* meExtraPhiAtBTLEff_;
  MonitorElement* meExtraMTDfailExtenderEtaEff_;
  MonitorElement* meExtraMTDfailExtenderPtEff_;
};

// ------------ constructor and destructor --------------
MtdTracksHarvester::MtdTracksHarvester(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")) {}

MtdTracksHarvester::~MtdTracksHarvester() {}

// auxiliary method to compute efficiency from the ratio of two 1D MonitorElement
void MtdTracksHarvester::computeEfficiency1D(MonitorElement* num, MonitorElement* den, MonitorElement* result) {
  for (int ibin = 1; ibin <= den->getNbinsX(); ibin++) {
    double eff = num->getBinContent(ibin) / den->getBinContent(ibin);
    double bin_err = sqrt((num->getBinContent(ibin) * (den->getBinContent(ibin) - num->getBinContent(ibin))) /
                          pow(den->getBinContent(ibin), 3));
    if (den->getBinContent(ibin) == 0) {
      eff = 0;
      bin_err = 0;
    }
    result->setBinContent(ibin, eff);
    result->setBinError(ibin, bin_err);
  }
}

void MtdTracksHarvester::normalize(MonitorElement* h, double scale) {
  double integral = h->getTH1F()->Integral();
  double norma = (integral > 0.) ? scale / integral : 0.;
  for (int ibin = 1; ibin <= h->getNbinsX(); ibin++) {
    double eff = h->getBinContent(ibin) * norma;
    double bin_err = h->getBinError(ibin) * norma;
    h->setBinContent(ibin, eff);
    h->setBinError(ibin, bin_err);
  }
}

// ------------ endjob tasks ----------------------------
void MtdTracksHarvester::dqmEndJob(DQMStore::IBooker& ibook, DQMStore::IGetter& igetter) {
  // --- Get the monitoring histograms
  MonitorElement* meBTLTrackEffEtaTot = igetter.get(folder_ + "TrackBTLEffEtaTot");
  MonitorElement* meBTLTrackEffPhiTot = igetter.get(folder_ + "TrackBTLEffPhiTot");
  MonitorElement* meBTLTrackEffPtTot = igetter.get(folder_ + "TrackBTLEffPtTot");
  MonitorElement* meBTLTrackEffEtaMtd = igetter.get(folder_ + "TrackBTLEffEtaMtd");
  MonitorElement* meBTLTrackEffPhiMtd = igetter.get(folder_ + "TrackBTLEffPhiMtd");
  MonitorElement* meBTLTrackEffPtMtd = igetter.get(folder_ + "TrackBTLEffPtMtd");

  MonitorElement* meETLTrackEffEtaTot = igetter.get(folder_ + "TrackETLEffEtaTot");
  MonitorElement* meETLTrackEffPhiTot = igetter.get(folder_ + "TrackETLEffPhiTot");
  MonitorElement* meETLTrackEffPtTot = igetter.get(folder_ + "TrackETLEffPtTot");
  MonitorElement* meETLTrackEffEtaMtd = igetter.get(folder_ + "TrackETLEffEtaMtd");
  MonitorElement* meETLTrackEffPhiMtd = igetter.get(folder_ + "TrackETLEffPhiMtd");
  MonitorElement* meETLTrackEffPtMtd = igetter.get(folder_ + "TrackETLEffPtMtd");
  MonitorElement* meETLTrackEffEta2Mtd = igetter.get(folder_ + "TrackETLEffEta2Mtd");
  MonitorElement* meETLTrackEffPhi2Mtd = igetter.get(folder_ + "TrackETLEffPhi2Mtd");
  MonitorElement* meETLTrackEffPt2Mtd = igetter.get(folder_ + "TrackETLEffPt2Mtd");

  MonitorElement* meETLTrackEffEtaTotLowPt0 = igetter.get(folder_ + "TrackETLEffEtaTotLowPt0");
  MonitorElement* meETLTrackEffEtaTotLowPt1 = igetter.get(folder_ + "TrackETLEffEtaTotLowPt1");
  MonitorElement* meETLTrackEffEtaMtdLowPt0 = igetter.get(folder_ + "TrackETLEffEtaMtdLowPt0");
  MonitorElement* meETLTrackEffEtaMtdLowPt1 = igetter.get(folder_ + "TrackETLEffEtaMtdLowPt1");
  MonitorElement* meETLTrackEffEta2MtdLowPt0 = igetter.get(folder_ + "TrackETLEffEta2MtdLowPt0");
  MonitorElement* meETLTrackEffEta2MtdLowPt1 = igetter.get(folder_ + "TrackETLEffEta2MtdLowPt1");
  
  MonitorElement* meTrackPtTot = igetter.get(folder_ + "TrackPtTot");
  MonitorElement* meExtraPtMtd = igetter.get(folder_ + "ExtraPtMtd");
  MonitorElement* meExtraPtEtl2Mtd = igetter.get(folder_ + "ExtraPtEtl2Mtd");
  MonitorElement* meTrackMatchedTPEffPtTot = igetter.get(folder_ + "MatchedTPEffPtTot");
  MonitorElement* meTrackMatchedTPEffPtTotLV = igetter.get(folder_ + "MatchedTPEffPtTotLV");
  MonitorElement* meTrackMatchedTPEffPtMtd = igetter.get(folder_ + "MatchedTPEffPtMtd");
  MonitorElement* meTrackMatchedTPEffPtEtl2Mtd = igetter.get(folder_ + "MatchedTPEffPtEtl2Mtd");
  MonitorElement* meTrackEtaTot = igetter.get(folder_ + "TrackEtaTot");
  MonitorElement* meExtraEtaMtd = igetter.get(folder_ + "ExtraEtaMtd");
  MonitorElement* meExtraEtaEtl2Mtd = igetter.get(folder_ + "ExtraEtaEtl2Mtd");
  MonitorElement* meTrackMatchedTPEffEtaTot = igetter.get(folder_ + "MatchedTPEffEtaTot");
  MonitorElement* meTrackMatchedTPEffEtaTotLV = igetter.get(folder_ + "MatchedTPEffEtaTotLV");
  MonitorElement* meTrackMatchedTPEffEtaMtd = igetter.get(folder_ + "MatchedTPEffEtaMtd");
  MonitorElement* meTrackMatchedTPEffEtaEtl2Mtd = igetter.get(folder_ + "MatchedTPEffEtaEtl2Mtd");

  //
  MonitorElement* meBTLTrackMatchedTPmtdDirectEta = igetter.get(folder_ + "BTLTrackMatchedTPmtdDirectEta");
  MonitorElement* meBTLTrackMatchedTPmtdDirectPt = igetter.get(folder_ + "BTLTrackMatchedTPmtdDirectPt");
  MonitorElement* meBTLTrackMatchedTPmtdOtherEta = igetter.get(folder_ + "BTLTrackMatchedTPmtdOtherEta"); 
  MonitorElement* meBTLTrackMatchedTPmtdOtherPt = igetter.get(folder_ + "BTLTrackMatchedTPmtdOtherPt"); ; 
  MonitorElement* meBTLTrackMatchedTPnomtdEta = igetter.get(folder_ + "BTLTrackMatchedTPnomtdEta");
  MonitorElement* meBTLTrackMatchedTPnomtdPt = igetter.get(folder_ + "BTLTrackMatchedTPnomtdPt");

  MonitorElement* meBTLTrackMatchedTPmtdDirectCorrectAssocEta = igetter.get(folder_ + "BTLTrackMatchedTPmtdDirectCorrectAssocEta");
  MonitorElement* meBTLTrackMatchedTPmtdDirectCorrectAssocPt = igetter.get(folder_ + "BTLTrackMatchedTPmtdDirectCorrectAssocPt");
  MonitorElement* meBTLTrackMatchedTPmtdDirectWrongAssocEta = igetter.get(folder_ + "BTLTrackMatchedTPmtdDirectWrongAssocEta");
  MonitorElement* meBTLTrackMatchedTPmtdDirectWrongAssocPt = igetter.get(folder_ + "BTLTrackMatchedTPmtdDirectWrongAssocPt");
  MonitorElement* meBTLTrackMatchedTPmtdDirectNoAssocEta = igetter.get(folder_ + "BTLTrackMatchedTPmtdDirectNoAssocEta");
  MonitorElement* meBTLTrackMatchedTPmtdDirectNoAssocPt = igetter.get(folder_ + "BTLTrackMatchedTPmtdDirectNoAssocPt");

  MonitorElement* meBTLTrackMatchedTPmtdOtherCorrectAssocEta = igetter.get(folder_ + "BTLTrackMatchedTPmtdOtherCorrectAssocEta");
  MonitorElement* meBTLTrackMatchedTPmtdOtherCorrectAssocPt = igetter.get(folder_ + "BTLTrackMatchedTPmtdOtherCorrectAssocPt");
  MonitorElement* meBTLTrackMatchedTPmtdOtherWrongAssocEta = igetter.get(folder_ + "BTLTrackMatchedTPmtdOtherWrongAssocEta");
  MonitorElement* meBTLTrackMatchedTPmtdOtherWrongAssocPt = igetter.get(folder_ + "BTLTrackMatchedTPmtdOtherWrongAssocPt");
  MonitorElement* meBTLTrackMatchedTPmtdOtherNoAssocEta = igetter.get(folder_ + "BTLTrackMatchedTPmtdOtherNoAssocEta");
  MonitorElement* meBTLTrackMatchedTPmtdOtherNoAssocPt = igetter.get(folder_ + "BTLTrackMatchedTPmtdOtherNoAssocPt");
  
  MonitorElement* meBTLTrackMatchedTPnomtdAssocEta = igetter.get(folder_ + "BTLTrackMatchedTPnomtdAssocEta");
  MonitorElement* meBTLTrackMatchedTPnomtdAssocPt = igetter.get(folder_ + "BTLTrackMatchedTPnomtdAssocPt");

  MonitorElement* meETLTrackMatchedTPmtd1Eta = igetter.get(folder_ + "ETLTrackMatchedTPmtd1Eta");
  MonitorElement* meETLTrackMatchedTPmtd1Pt = igetter.get(folder_ + "ETLTrackMatchedTPmtd1Pt");
  MonitorElement* meETLTrackMatchedTPmtd2Eta = igetter.get(folder_ + "ETLTrackMatchedTPmtd2Eta"); 
  MonitorElement* meETLTrackMatchedTPmtd2Pt = igetter.get(folder_ + "ETLTrackMatchedTPmtd2Pt"); ; 
  MonitorElement* meETLTrackMatchedTPnomtdEta = igetter.get(folder_ + "ETLTrackMatchedTPnomtdEta");
  MonitorElement* meETLTrackMatchedTPnomtdPt = igetter.get(folder_ + "ETLTrackMatchedTPnomtdPt");

  MonitorElement* meETLTrackMatchedTPmtd1CorrectAssocEta = igetter.get(folder_ + "ETLTrackMatchedTPmtd1CorrectAssocEta");
  MonitorElement* meETLTrackMatchedTPmtd1CorrectAssocPt = igetter.get(folder_ + "ETLTrackMatchedTPmtd1CorrectAssocPt");
  MonitorElement* meETLTrackMatchedTPmtd1WrongAssocEta = igetter.get(folder_ + "ETLTrackMatchedTPmtd1WrongAssocEta");
  MonitorElement* meETLTrackMatchedTPmtd1WrongAssocPt = igetter.get(folder_ + "ETLTrackMatchedTPmtd1WrongAssocPt");
  MonitorElement* meETLTrackMatchedTPmtd1NoAssocEta = igetter.get(folder_ + "ETLTrackMatchedTPmtd1NoAssocEta");
  MonitorElement* meETLTrackMatchedTPmtd1NoAssocPt = igetter.get(folder_ + "ETLTrackMatchedTPmtd1NoAssocPt");

  MonitorElement* meETLTrackMatchedTPmtd2CorrectAssocEta = igetter.get(folder_ + "ETLTrackMatchedTPmtd2CorrectAssocEta");
  MonitorElement* meETLTrackMatchedTPmtd2CorrectAssocPt = igetter.get(folder_ + "ETLTrackMatchedTPmtd2CorrectAssocPt");
  MonitorElement* meETLTrackMatchedTPmtd2WrongAssocEta = igetter.get(folder_ + "ETLTrackMatchedTPmtd2WrongAssocEta");
  MonitorElement* meETLTrackMatchedTPmtd2WrongAssocPt = igetter.get(folder_ + "ETLTrackMatchedTPmtd2WrongAssocPt");
  MonitorElement* meETLTrackMatchedTPmtd2NoAssocEta = igetter.get(folder_ + "ETLTrackMatchedTPmtd2NoAssocEta");
  MonitorElement* meETLTrackMatchedTPmtd2NoAssocPt = igetter.get(folder_ + "ETLTrackMatchedTPmtd2NoAssocPt");
  
  MonitorElement* meETLTrackMatchedTPnomtdAssocEta = igetter.get(folder_ + "ETLTrackMatchedTPnomtdAssocEta");
  MonitorElement* meETLTrackMatchedTPnomtdAssocPt = igetter.get(folder_ + "ETLTrackMatchedTPnomtdAssocPt");
  
  //  
  MonitorElement* meTrackNumHits = igetter.get(folder_ + "TrackNumHits");
  MonitorElement* meTrackNumHitsNT = igetter.get(folder_ + "TrackNumHitsNT");
  MonitorElement* meExtraPhiAtBTL = igetter.get(folder_ + "ExtraPhiAtBTL");
  MonitorElement* meExtraPhiAtBTLmatched = igetter.get(folder_ + "ExtraPhiAtBTLmatched");
  MonitorElement* meExtraBTLeneInCone = igetter.get(folder_ + "ExtraBTLeneInCone");
  MonitorElement* meExtraMTDfailExtenderEta = igetter.get(folder_ + "ExtraMTDfailExtenderEta");
  MonitorElement* meExtraMTDfailExtenderPt = igetter.get(folder_ + "ExtraMTDfailExtenderPt");

  if (!meBTLTrackEffEtaTot || !meBTLTrackEffPhiTot || !meBTLTrackEffPtTot ||
      !meBTLTrackEffEtaMtd || !meBTLTrackEffPhiMtd || !meBTLTrackEffPtMtd ||
      !meETLTrackEffEtaTot || !meETLTrackEffPhiTot || !meETLTrackEffPtTot ||
      !meETLTrackEffEtaMtd || !meETLTrackEffPhiMtd || !meETLTrackEffPtMtd ||
      !meETLTrackEffEta2Mtd || !meETLTrackEffPhi2Mtd || !meETLTrackEffPt2Mtd ||
      !meETLTrackEffEtaTotLowPt0 || !meETLTrackEffEtaTotLowPt1 ||
      !meETLTrackEffEtaMtdLowPt0 || !meETLTrackEffEtaMtdLowPt1 || 
      !meETLTrackEffEta2MtdLowPt0 || !meETLTrackEffEta2MtdLowPt1 ||
      !meTrackMatchedTPEffPtTotLV || !meTrackMatchedTPEffPtMtd || !meTrackMatchedTPEffPtEtl2Mtd ||
      !meTrackMatchedTPEffEtaTot ||
      !meTrackMatchedTPEffEtaTotLV ||!meTrackMatchedTPEffEtaMtd || !meTrackMatchedTPEffEtaEtl2Mtd ||
      !meBTLTrackMatchedTPmtdDirectEta || !meBTLTrackMatchedTPmtdDirectPt || !meBTLTrackMatchedTPmtdOtherEta || !meBTLTrackMatchedTPmtdOtherPt || !meBTLTrackMatchedTPnomtdEta || !meBTLTrackMatchedTPnomtdPt ||
      !meBTLTrackMatchedTPmtdDirectCorrectAssocEta || !meBTLTrackMatchedTPmtdDirectCorrectAssocPt || !meBTLTrackMatchedTPmtdDirectWrongAssocEta || !meBTLTrackMatchedTPmtdDirectWrongAssocPt || !meBTLTrackMatchedTPmtdDirectNoAssocEta || !meBTLTrackMatchedTPmtdDirectNoAssocPt ||
      !meBTLTrackMatchedTPmtdOtherCorrectAssocEta || !meBTLTrackMatchedTPmtdOtherCorrectAssocPt || !meBTLTrackMatchedTPmtdOtherWrongAssocEta || !meBTLTrackMatchedTPmtdOtherWrongAssocPt || !meBTLTrackMatchedTPmtdOtherNoAssocEta || !meBTLTrackMatchedTPmtdOtherNoAssocPt ||
      !meBTLTrackMatchedTPnomtdAssocEta || !meBTLTrackMatchedTPnomtdAssocPt ||
      !meETLTrackMatchedTPmtd1Eta || !meETLTrackMatchedTPmtd1Pt || !meETLTrackMatchedTPmtd2Eta || !meETLTrackMatchedTPmtd2Pt || !meETLTrackMatchedTPnomtdEta || !meETLTrackMatchedTPnomtdPt ||
      !meETLTrackMatchedTPmtd1CorrectAssocEta || !meETLTrackMatchedTPmtd1CorrectAssocPt || !meETLTrackMatchedTPmtd1WrongAssocEta || !meETLTrackMatchedTPmtd1WrongAssocPt || !meETLTrackMatchedTPmtd1NoAssocEta || !meETLTrackMatchedTPmtd1NoAssocPt ||
      !meETLTrackMatchedTPmtd2CorrectAssocEta || !meETLTrackMatchedTPmtd2CorrectAssocPt || !meETLTrackMatchedTPmtd2WrongAssocEta || !meETLTrackMatchedTPmtd2WrongAssocPt || !meETLTrackMatchedTPmtd2NoAssocEta || !meETLTrackMatchedTPmtd2NoAssocPt ||
      !meETLTrackMatchedTPnomtdAssocEta || !meETLTrackMatchedTPnomtdAssocPt ||
      !meTrackNumHits || !meTrackNumHitsNT ||
      !meTrackPtTot || !meTrackEtaTot || !meExtraPtMtd || !meExtraPtEtl2Mtd || !meExtraEtaMtd || !meExtraEtaEtl2Mtd ||
      !meExtraPhiAtBTL || !meExtraPhiAtBTLmatched || !meExtraBTLeneInCone || !meExtraMTDfailExtenderEta ||
      !meExtraMTDfailExtenderPt) {
    edm::LogError("MtdTracksHarvester") << "Monitoring histograms not found!" << std::endl;
    return;
  }

  // --- Book  histograms
  ibook.cd(folder_);
  meBtlEtaEff_ = ibook.book1D("BtlEtaEff",
                              " Track Efficiency VS Eta;#eta;Efficiency",
                              meBTLTrackEffEtaTot->getNbinsX(),
                              meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                              meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBtlEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackEffEtaMtd, meBTLTrackEffEtaTot, meBtlEtaEff_);

  meBtlPhiEff_ = ibook.book1D("BtlPhiEff",
                              "Track Efficiency VS Phi;#phi [rad];Efficiency",
                              meBTLTrackEffPhiTot->getNbinsX(),
                              meBTLTrackEffPhiTot->getTH1()->GetXaxis()->GetXmin(),
                              meBTLTrackEffPhiTot->getTH1()->GetXaxis()->GetXmax());
  meBtlPhiEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackEffPhiMtd, meBTLTrackEffPhiTot, meBtlPhiEff_);

  meBtlPtEff_ = ibook.book1D("BtlPtEff",
                             "Track Efficiency VS Pt;Pt [GeV];Efficiency",
                             meBTLTrackEffPtTot->getNbinsX(),
                             meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                             meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBtlPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackEffPtMtd, meBTLTrackEffPtTot, meBtlPtEff_);

  meEtlEtaEff_ = ibook.book1D("EtlEtaEff",
			      " Track Efficiency VS Eta;#eta;Efficiency",
			      meETLTrackEffEtaTot->getNbinsX(),
			      meETLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
			      meETLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEtaMtd, meETLTrackEffEtaTot, meEtlEtaEff_);

  meEtlPhiEff_ = ibook.book1D("EtlPhiEff",
			      "Track Efficiency VS Phi;#phi [rad];Efficiency",
			      meETLTrackEffPhiTot->getNbinsX(),
			      meETLTrackEffPhiTot->getTH1()->GetXaxis()->GetXmin(),
			      meETLTrackEffPhiTot->getTH1()->GetXaxis()->GetXmax());
  meEtlPhiEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPhiMtd, meETLTrackEffPhiTot, meEtlPhiEff_);

  meEtlPtEff_ = ibook.book1D("EtlPtEff",
			     "Track Efficiency VS Pt;Pt [GeV];Efficiency",
			     meETLTrackEffPtTot->getNbinsX(),
			     meETLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
			     meETLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meEtlPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPtMtd, meETLTrackEffPtTot, meEtlPtEff_);
  
  meEtlEtaEff2_ = ibook.book1D("EtlEtaEff2",
			       " Track Efficiency VS Eta (2 hits);#eta;Efficiency",
			       meETLTrackEffEtaTot->getNbinsX(),
			       meETLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
			       meETLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEff2_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEta2Mtd, meETLTrackEffEtaTot, meEtlEtaEff2_);
  
  meEtlPhiEff2_ = ibook.book1D("EtlPhiEff2",
			       "Track Efficiency VS Phi (2 hits);#phi [rad];Efficiency",
                                  meETLTrackEffPhiTot->getNbinsX(),
                                  meETLTrackEffPhiTot->getTH1()->GetXaxis()->GetXmin(),
                                  meETLTrackEffPhiTot->getTH1()->GetXaxis()->GetXmax());
  meEtlPhiEff2_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPhi2Mtd, meETLTrackEffPhiTot, meEtlPhiEff2_);

  meEtlPtEff2_ = ibook.book1D("EtlPtEff2",
			      "Track Efficiency VS Pt (2 hits);Pt [GeV];Efficiency",
			      meETLTrackEffPtTot->getNbinsX(),
			      meETLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
			      meETLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meEtlPtEff2_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPt2Mtd, meETLTrackEffPtTot, meEtlPtEff2_);

  // low pT
  meEtlEtaEffLowPt_[0] = ibook.book1D("EtlEtaEffLowPt0",
                                      " Track Efficiency VS Eta, 0.2 < pt < 0.45;#eta;Efficiency",
                                      meETLTrackEffEtaTotLowPt0->getNbinsX(),
                                      meETLTrackEffEtaTotLowPt0->getTH1()->GetXaxis()->GetXmin(),
                                      meETLTrackEffEtaTotLowPt0->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEffLowPt_[0]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEtaMtdLowPt0, meETLTrackEffEtaTotLowPt0, meEtlEtaEffLowPt_[0]);
  
  meEtlEtaEffLowPt_[1] = ibook.book1D("EtlEtaEffLowPt1",
                                      " Track Efficiency VS Eta, 0.45 < pt < 0.7;#eta;Efficiency",
                                      meETLTrackEffEtaTotLowPt1->getNbinsX(),
                                      meETLTrackEffEtaTotLowPt1->getTH1()->GetXaxis()->GetXmin(),
                                      meETLTrackEffEtaTotLowPt1->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEffLowPt_[1]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEtaMtdLowPt1, meETLTrackEffEtaTotLowPt1, meEtlEtaEffLowPt_[1]);

  meEtlEtaEff2LowPt_[0] = ibook.book1D("EtlEtaEff2LowPt0",
                                       " Track Efficiency VS Eta (2 hits), 0.2 < pt < 0.45;#eta;Efficiency",
                                       meETLTrackEffEtaTotLowPt0->getNbinsX(),
                                       meETLTrackEffEtaTotLowPt0->getTH1()->GetXaxis()->GetXmin(),
                                       meETLTrackEffEtaTotLowPt0->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEff2LowPt_[0]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEta2MtdLowPt0, meETLTrackEffEtaTotLowPt0, meEtlEtaEff2LowPt_[0]);


  meEtlEtaEff2LowPt_[1] = ibook.book1D("EtlEtaEff2LowPt1",
                                       " Track Efficiency VS Eta (2 hits), 0.45 < pt < 0.7;#eta;Efficiency",
                                       meETLTrackEffEtaTotLowPt1->getNbinsX(),
                                       meETLTrackEffEtaTotLowPt1->getTH1()->GetXaxis()->GetXmin(),
                                       meETLTrackEffEtaTotLowPt1->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEff2LowPt_[1]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEta2MtdLowPt1, meETLTrackEffEtaTotLowPt1, meEtlEtaEff2LowPt_[1]);
  

  meExtraPtEff_ =
      ibook.book1D("ExtraPtEff",
                   "MTD matching efficiency wrt extrapolated track associated to LV VS Pt;Pt [GeV];Efficiency",
                   meTrackMatchedTPEffPtTotLV->getNbinsX(),
                   meTrackMatchedTPEffPtTotLV->getTH1()->GetXaxis()->GetXmin(),
                   meTrackMatchedTPEffPtTotLV->getTH1()->GetXaxis()->GetXmax());
  meExtraPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraPtMtd, meTrackMatchedTPEffPtTotLV, meExtraPtEff_);

  meExtraPtEtl2Eff_ =
      ibook.book1D("ExtraPtEtl2Eff",
                   "MTD matching efficiency (2 ETL) wrt extrapolated track associated to LV VS Pt;Pt [GeV];Efficiency",
                   meTrackMatchedTPEffPtTotLV->getNbinsX(),
                   meTrackMatchedTPEffPtTotLV->getTH1()->GetXaxis()->GetXmin(),
                   meTrackMatchedTPEffPtTotLV->getTH1()->GetXaxis()->GetXmax());
  meExtraPtEtl2Eff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraPtEtl2Mtd, meTrackMatchedTPEffPtTotLV, meExtraPtEtl2Eff_);

  meExtraEtaEff_ = ibook.book1D("ExtraEtaEff",
                                "MTD matching efficiency wrt extrapolated track associated to LV VS Eta;Eta;Efficiency",
                                meTrackMatchedTPEffEtaTotLV->getNbinsX(),
                                meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmin(),
                                meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmax());
  meExtraEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraEtaMtd, meTrackMatchedTPEffEtaTotLV, meExtraEtaEff_);

  meExtraEtaEtl2Eff_ =
      ibook.book1D("ExtraEtaEtl2Eff",
                   "MTD matching efficiency (2 ETL) wrt extrapolated track associated to LV VS Eta;Eta;Efficiency",
                   meTrackMatchedTPEffEtaTotLV->getNbinsX(),
                   meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmin(),
                   meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmax());
  meExtraEtaEtl2Eff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraEtaEtl2Mtd, meTrackMatchedTPEffEtaTotLV, meExtraEtaEtl2Eff_);

  meTPPtSelEff_ = ibook.book1D("TPPtSelEff",
                               "Track selected efficiency TP VS Pt;Pt [GeV];Efficiency",
                               meTrackPtTot->getNbinsX(),
                               meTrackPtTot->getTH1()->GetXaxis()->GetXmin(),
                               meTrackPtTot->getTH1()->GetXaxis()->GetXmax());
  meTPPtSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPEffPtTot, meTrackPtTot, meTPPtSelEff_);

  meTPEtaSelEff_ = ibook.book1D("TPEtaSelEff",
                                "Track selected efficiency TP VS Eta;Eta;Efficiency",
                                meTrackEtaTot->getNbinsX(),
                                meTrackEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                meTrackEtaTot->getTH1()->GetXaxis()->GetXmax());
  meTPEtaSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPEffEtaTot, meTrackEtaTot, meTPEtaSelEff_);

  meTPPtMatchEff_ = ibook.book1D("TPPtMatchEff",
                                 "Track matched to TP efficiency VS Pt;Pt [GeV];Efficiency",
                                 meTrackMatchedTPEffPtTot->getNbinsX(),
                                 meTrackMatchedTPEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                 meTrackMatchedTPEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meTPPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPEffPtMtd, meTrackMatchedTPEffPtTot, meTPPtMatchEff_);

  meTPEtaMatchEff_ = ibook.book1D("TPEtaMatchEff",
                                  "Track matched to TP efficiency VS Eta;Eta;Efficiency",
                                  meTrackMatchedTPEffEtaTot->getNbinsX(),
                                  meTrackMatchedTPEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                  meTrackMatchedTPEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meTPEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPEffEtaMtd, meTrackMatchedTPEffEtaTot, meTPEtaMatchEff_);

  meTPPtMatchEtl2Eff_ = ibook.book1D("TPPtMatchEtl2Eff",
                                     "Track matched to TP efficiency VS Pt, 2 ETL hits;Pt [GeV];Efficiency",
                                     meTrackMatchedTPEffPtTot->getNbinsX(),
                                     meTrackMatchedTPEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                     meTrackMatchedTPEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meTPPtMatchEtl2Eff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPEffPtEtl2Mtd, meTrackMatchedTPEffPtTot, meTPPtMatchEtl2Eff_);

  meTPEtaMatchEtl2Eff_ = ibook.book1D("TPEtaMatchEtl2Eff",
                                      "Track matched to TP efficiency VS Eta, 2 ETL hits;Eta;Efficiency",
                                      meTrackMatchedTPEffEtaTot->getNbinsX(),
                                      meTrackMatchedTPEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                      meTrackMatchedTPEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meTPEtaMatchEtl2Eff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPEffEtaEtl2Mtd, meTrackMatchedTPEffEtaTot, meTPEtaMatchEtl2Eff_);


  // == Track-cluster matching efficiencies based on mc truth
  // -- BTL
  meBTLTPmtdDirectEtaSelEff_ = ibook.book1D("BTLTPmtdDirectEtaSelEff",
					    "Track selected efficiency TP-mtd hit (direct) VS Eta",
					    meBTLTrackEffEtaTot->getNbinsX(),
					    meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
					    meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdDirectEtaSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdDirectEta, meBTLTrackEffEtaTot, meBTLTPmtdDirectEtaSelEff_);
  
  meBTLTPmtdDirectPtSelEff_ = ibook.book1D("BTLTPmtdDirectPtSelEff",
					   "Track selected efficiency TP-mtd hit (direct) VS Pt",
					   meBTLTrackEffPtTot->getNbinsX(),
					   meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
					   meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdDirectPtSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdDirectPt, meBTLTrackEffPtTot, meBTLTPmtdDirectPtSelEff_);
  
  meBTLTPmtdOtherEtaSelEff_ = ibook.book1D("BTLTPmtdOtherEtaSelEff",
					   "Track selected efficiency TP-mtd hit (other) VS Eta",
					   meBTLTrackEffEtaTot->getNbinsX(),
					   meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
					   meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdOtherEtaSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdOtherEta, meBTLTrackEffEtaTot, meBTLTPmtdOtherEtaSelEff_);
  
  meBTLTPmtdOtherPtSelEff_ = ibook.book1D("BTLTPmtdOtherPtSelEff",
					  "Track selected efficiency TP-mtd hit (other) VS Pt",
					  meBTLTrackEffPtTot->getNbinsX(),
					  meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
					  meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdOtherPtSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdOtherPt, meBTLTrackEffPtTot, meBTLTPmtdOtherPtSelEff_);		     
  
  meBTLTPnomtdEtaSelEff_ = ibook.book1D("BTLTPnomtdEtaSelEff",
					"Track selected efficiency TP-no mtd hit VS Eta",
					meBTLTrackEffEtaTot->getNbinsX(),
					meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
					meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLTPnomtdEtaSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPnomtdEta, meBTLTrackEffEtaTot, meBTLTPnomtdEtaSelEff_);

  meBTLTPnomtdPtSelEff_ = ibook.book1D("BTLTPnomtdPtSelEff",
				      "Track selected efficiency TP-no mtd hit VS Pt",
				      meBTLTrackEffPtTot->getNbinsX(),
				      meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
				      meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLTPnomtdPtSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPnomtdPt, meBTLTrackEffPtTot, meBTLTPnomtdPtSelEff_);

    
  meBTLTPmtdDirectCorrectAssocEtaMatchEff_ = ibook.book1D("BTLTPmtdDirectCorrectAssocEtaMatchEff",
						       "Track efficiency TP-mtd hit (direct), correct reco match VS Eta",
						       meBTLTrackMatchedTPmtdDirectEta->getNbinsX(),
						       meBTLTrackMatchedTPmtdDirectEta->getTH1()->GetXaxis()->GetXmin(),
						       meBTLTrackMatchedTPmtdDirectEta->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdDirectCorrectAssocEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdDirectCorrectAssocEta, meBTLTrackMatchedTPmtdDirectEta, meBTLTPmtdDirectCorrectAssocEtaMatchEff_);

  meBTLTPmtdDirectCorrectAssocPtMatchEff_ = ibook.book1D("BTLTPmtdDirectCorrectAssocPtMatchEff",
						      "Track efficiency TP-mtd hit (direct), correct reco match VS Pt",
						      meBTLTrackMatchedTPmtdDirectPt->getNbinsX(),
						      meBTLTrackMatchedTPmtdDirectPt->getTH1()->GetXaxis()->GetXmin(),
						      meBTLTrackMatchedTPmtdDirectPt->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdDirectCorrectAssocPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdDirectCorrectAssocPt, meBTLTrackMatchedTPmtdDirectPt, meBTLTPmtdDirectCorrectAssocPtMatchEff_);
  
  meBTLTPmtdDirectWrongAssocEtaMatchEff_ = ibook.book1D("BTLTPmtdDirectWrongAssocEtaMatchEff",
						     "Track efficiency TP-mtd hit (direct), incorrect reco match VS Eta",
						     meBTLTrackMatchedTPmtdDirectEta->getNbinsX(),
						     meBTLTrackMatchedTPmtdDirectEta->getTH1()->GetXaxis()->GetXmin(),
						     meBTLTrackMatchedTPmtdDirectEta->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdDirectWrongAssocEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdDirectWrongAssocEta, meBTLTrackMatchedTPmtdDirectEta, meBTLTPmtdDirectWrongAssocEtaMatchEff_);
  
  meBTLTPmtdDirectWrongAssocPtMatchEff_ = ibook.book1D("BTLTPmtdDirectWrongAssocPtMatchEff",
						    "Track efficiency TP-mtd hit (direct), incorrect reco match VS Pt",
						    meBTLTrackMatchedTPmtdDirectPt->getNbinsX(),
						    meBTLTrackMatchedTPmtdDirectPt->getTH1()->GetXaxis()->GetXmin(),
						    meBTLTrackMatchedTPmtdDirectPt->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdDirectWrongAssocPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdDirectWrongAssocPt, meBTLTrackMatchedTPmtdDirectPt, meBTLTPmtdDirectWrongAssocPtMatchEff_);

  
  meBTLTPmtdDirectNoAssocEtaMatchEff_ = ibook.book1D("BTLTPmtdDirectNoAssocEtaMatchEff",
						     "Track efficiency TP-mtd hit (direct), no reco match VS Eta",
						     meBTLTrackMatchedTPmtdDirectEta->getNbinsX(),
						     meBTLTrackMatchedTPmtdDirectEta->getTH1()->GetXaxis()->GetXmin(),
						     meBTLTrackMatchedTPmtdDirectEta->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdDirectNoAssocEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdDirectNoAssocEta, meBTLTrackMatchedTPmtdDirectEta, meBTLTPmtdDirectNoAssocEtaMatchEff_);

  meBTLTPmtdDirectNoAssocPtMatchEff_ = ibook.book1D("BTLTPmtdDirectNoAssocPtMatchEff",
						    "Track efficiency TP-mtd hit (direct), no reco match VS Pt",
						    meBTLTrackMatchedTPmtdDirectPt->getNbinsX(),
						    meBTLTrackMatchedTPmtdDirectPt->getTH1()->GetXaxis()->GetXmin(),
						    meBTLTrackMatchedTPmtdDirectPt->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdDirectNoAssocPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdDirectNoAssocPt, meBTLTrackMatchedTPmtdDirectPt, meBTLTPmtdDirectNoAssocPtMatchEff_);


  
  meBTLTPmtdOtherCorrectAssocEtaMatchEff_ = ibook.book1D("BTLTPmtdOtherCorrectAssocEtaMatchEff",
						       "Track efficiency TP-mtd hit (other), correct reco match VS Eta",
						       meBTLTrackMatchedTPmtdOtherEta->getNbinsX(),
						       meBTLTrackMatchedTPmtdOtherEta->getTH1()->GetXaxis()->GetXmin(),
						       meBTLTrackMatchedTPmtdOtherEta->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdOtherCorrectAssocEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdOtherCorrectAssocEta, meBTLTrackMatchedTPmtdOtherEta, meBTLTPmtdOtherCorrectAssocEtaMatchEff_);

  meBTLTPmtdOtherCorrectAssocPtMatchEff_ = ibook.book1D("BTLTPmtdOtherCorrectAssocPtMatchEff",
						      "Track efficiency TP-mtd hit (other), correct reco match VS Pt",
						      meBTLTrackMatchedTPmtdOtherPt->getNbinsX(),
						      meBTLTrackMatchedTPmtdOtherPt->getTH1()->GetXaxis()->GetXmin(),
						      meBTLTrackMatchedTPmtdOtherPt->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdOtherCorrectAssocPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdOtherCorrectAssocPt, meBTLTrackMatchedTPmtdOtherPt, meBTLTPmtdOtherCorrectAssocPtMatchEff_);
  
  meBTLTPmtdOtherWrongAssocEtaMatchEff_ = ibook.book1D("BTLTPmtdOtherWrongAssocEtaMatchEff",
						     "Track efficiency TP-mtd hit (other), incorrect reco match VS Eta",
						     meBTLTrackMatchedTPmtdOtherEta->getNbinsX(),
						     meBTLTrackMatchedTPmtdOtherEta->getTH1()->GetXaxis()->GetXmin(),
						     meBTLTrackMatchedTPmtdOtherEta->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdOtherWrongAssocEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdOtherWrongAssocEta, meBTLTrackMatchedTPmtdOtherEta, meBTLTPmtdOtherWrongAssocEtaMatchEff_);
  
  meBTLTPmtdOtherWrongAssocPtMatchEff_ = ibook.book1D("BTLTPmtdOtherWrongAssocPtMatchEff",
						    "Track efficiency TP-mtd hit (other), incorrect reco match VS Pt",
						    meBTLTrackMatchedTPmtdOtherPt->getNbinsX(),
						    meBTLTrackMatchedTPmtdOtherPt->getTH1()->GetXaxis()->GetXmin(),
						    meBTLTrackMatchedTPmtdOtherPt->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdOtherWrongAssocPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdOtherWrongAssocPt, meBTLTrackMatchedTPmtdOtherPt, meBTLTPmtdOtherWrongAssocPtMatchEff_);

  
  meBTLTPmtdOtherNoAssocEtaMatchEff_ = ibook.book1D("BTLTPmtdOtherNoAssocEtaMatchEff",
						 "Track efficiency TP-mtd hit (other), no reco match VS Eta",
						 meBTLTrackMatchedTPmtdOtherEta->getNbinsX(),
						 meBTLTrackMatchedTPmtdOtherEta->getTH1()->GetXaxis()->GetXmin(),
						 meBTLTrackMatchedTPmtdOtherEta->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdOtherNoAssocEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdOtherNoAssocEta, meBTLTrackMatchedTPmtdOtherEta, meBTLTPmtdOtherNoAssocEtaMatchEff_);
  
  meBTLTPmtdOtherNoAssocPtMatchEff_ = ibook.book1D("BTLTPmtdOtherNoAssocPtMatchEff",
						"Track efficiency TP-mtd hit (other), no reco match VS Pt",
						meBTLTrackMatchedTPmtdOtherPt->getNbinsX(),
						meBTLTrackMatchedTPmtdOtherPt->getTH1()->GetXaxis()->GetXmin(),
						meBTLTrackMatchedTPmtdOtherPt->getTH1()->GetXaxis()->GetXmax());
  meBTLTPmtdOtherNoAssocPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPmtdOtherNoAssocPt, meBTLTrackMatchedTPmtdOtherPt, meBTLTPmtdOtherNoAssocPtMatchEff_);

  
  meBTLTPnomtdEtaMatchEff_ = ibook.book1D("BTLTPnomtdEtaMatchEff",
					    "Track efficiency TP- no mtd hit, with reco match VS Eta",
					    meBTLTrackMatchedTPnomtdEta->getNbinsX(),
					    meBTLTrackMatchedTPnomtdEta->getTH1()->GetXaxis()->GetXmin(),
					    meBTLTrackMatchedTPnomtdEta->getTH1()->GetXaxis()->GetXmax());
  meBTLTPnomtdEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPnomtdAssocEta, meBTLTrackMatchedTPnomtdEta, meBTLTPnomtdEtaMatchEff_);

  meBTLTPnomtdPtMatchEff_ = ibook.book1D("BTLTPnomtdPtMatchEff",
					    "Track efficiency TP- no mtd hit, with reco match VS Pt",
					    meBTLTrackMatchedTPnomtdPt->getNbinsX(),
					    meBTLTrackMatchedTPnomtdPt->getTH1()->GetXaxis()->GetXmin(),
					    meBTLTrackMatchedTPnomtdPt->getTH1()->GetXaxis()->GetXmax());
  meBTLTPnomtdPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackMatchedTPnomtdAssocPt, meBTLTrackMatchedTPnomtdPt, meBTLTPnomtdPtMatchEff_);


  // -- ETL
  meETLTPmtd1EtaSelEff_ = ibook.book1D("ETLTPmtd1EtaSelEff",
				       "Track selected efficiency TP-mtd hit (>=1 sim hit) VS Eta",
				       meETLTrackEffEtaTot->getNbinsX(),
				       meETLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
				       meETLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd1EtaSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd1Eta, meETLTrackEffEtaTot, meETLTPmtd1EtaSelEff_);
  
  meETLTPmtd1PtSelEff_ = ibook.book1D("ETLTPmtd1PtSelEff",
				      "Track selected efficiency TP-mtd hit (>=1 sim hit) VS Pt",
				      meETLTrackEffPtTot->getNbinsX(),
				      meETLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
				      meETLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd1PtSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd1Pt, meETLTrackEffPtTot, meETLTPmtd1PtSelEff_);
  
  meETLTPmtd2EtaSelEff_ = ibook.book1D("ETLTPmtd2EtaSelEff",
				       "Track selected efficiency TP-mtd hit (2 sim hits) VS Eta",
				       meETLTrackEffEtaTot->getNbinsX(),
				       meETLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
				       meETLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd2EtaSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd2Eta, meETLTrackEffEtaTot, meETLTPmtd2EtaSelEff_);
  
  meETLTPmtd2PtSelEff_ = ibook.book1D("ETLTPmtd2PtSelEff",
				      "Track selected efficiency TP-mtd hit (2 sim hits) VS Pt",
				      meETLTrackEffPtTot->getNbinsX(),
				      meETLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
				      meETLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd2PtSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd2Pt, meETLTrackEffPtTot, meETLTPmtd2PtSelEff_);		     
  
  meETLTPnomtdEtaSelEff_ = ibook.book1D("ETLTPnomtdEtaSelEff",
					"Track selected efficiency TP-no mtd hit VS Eta",
					meETLTrackEffEtaTot->getNbinsX(),
					meETLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
					meETLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meETLTPnomtdEtaSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPnomtdEta, meETLTrackEffEtaTot, meETLTPnomtdEtaSelEff_);
  
  meETLTPnomtdPtSelEff_ = ibook.book1D("ETLTPnomtdPtSelEff",
				       "Track selected efficiency TP-no mtd hit VS Pt",
				       meETLTrackEffPtTot->getNbinsX(),
				       meETLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
				       meETLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meETLTPnomtdPtSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPnomtdPt, meETLTrackEffPtTot, meETLTPnomtdPtSelEff_);
  
  
  meETLTPmtd1CorrectAssocEtaMatchEff_ = ibook.book1D("ETLTPmtd1CorrectAssocEtaMatchEff",
						     "Track efficiency TP-mtd hit (>=1 sim hit), correct reco match VS Eta",
						     meETLTrackMatchedTPmtd1Eta->getNbinsX(),
						     meETLTrackMatchedTPmtd1Eta->getTH1()->GetXaxis()->GetXmin(),
						     meETLTrackMatchedTPmtd1Eta->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd1CorrectAssocEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd1CorrectAssocEta, meETLTrackMatchedTPmtd1Eta, meETLTPmtd1CorrectAssocEtaMatchEff_);
  
  meETLTPmtd1CorrectAssocPtMatchEff_ = ibook.book1D("ETLTPmtd1CorrectAssocPtMatchEff",
						    "Track efficiency TP-mtd hit (>=1 sim hit), correct reco match VS Pt",
						    meETLTrackMatchedTPmtd1Pt->getNbinsX(),
						    meETLTrackMatchedTPmtd1Pt->getTH1()->GetXaxis()->GetXmin(),
						    meETLTrackMatchedTPmtd1Pt->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd1CorrectAssocPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd1CorrectAssocPt, meETLTrackMatchedTPmtd1Pt, meETLTPmtd1CorrectAssocPtMatchEff_);
  
  meETLTPmtd1WrongAssocEtaMatchEff_ = ibook.book1D("ETLTPmtd1WrongAssocEtaMatchEff",
						   "Track efficiency TP-mtd hit (>=1 sim hit), incorrect reco match VS Eta",
						   meETLTrackMatchedTPmtd1Eta->getNbinsX(),
						   meETLTrackMatchedTPmtd1Eta->getTH1()->GetXaxis()->GetXmin(),
						   meETLTrackMatchedTPmtd1Eta->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd1WrongAssocEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd1WrongAssocEta, meETLTrackMatchedTPmtd1Eta, meETLTPmtd1WrongAssocEtaMatchEff_);
  
  meETLTPmtd1WrongAssocPtMatchEff_ = ibook.book1D("ETLTPmtd1WrongAssocPtMatchEff",
						  "Track efficiency TP-mtd hit (>=1 sim hit), incorrect reco match VS Pt",
						  meETLTrackMatchedTPmtd1Pt->getNbinsX(),
						  meETLTrackMatchedTPmtd1Pt->getTH1()->GetXaxis()->GetXmin(),
						  meETLTrackMatchedTPmtd1Pt->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd1WrongAssocPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd1WrongAssocPt, meETLTrackMatchedTPmtd1Pt, meETLTPmtd1WrongAssocPtMatchEff_);
  
  
  meETLTPmtd1NoAssocEtaMatchEff_ = ibook.book1D("ETLTPmtd1NoAssocEtaMatchEff",
						"Track efficiency TP-mtd hit (>=1 sim hit), no reco match VS Eta",
						meETLTrackMatchedTPmtd1Eta->getNbinsX(),
						meETLTrackMatchedTPmtd1Eta->getTH1()->GetXaxis()->GetXmin(),
						meETLTrackMatchedTPmtd1Eta->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd1NoAssocEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd1NoAssocEta, meETLTrackMatchedTPmtd1Eta, meETLTPmtd1NoAssocEtaMatchEff_);
  
  meETLTPmtd1NoAssocPtMatchEff_ = ibook.book1D("ETLTPmtd1NoAssocPtMatchEff",
					       "Track efficiency TP-mtd hit (>=1 sim hit), no reco match VS Pt",
					       meETLTrackMatchedTPmtd1Pt->getNbinsX(),
					       meETLTrackMatchedTPmtd1Pt->getTH1()->GetXaxis()->GetXmin(),
					       meETLTrackMatchedTPmtd1Pt->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd1NoAssocPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd1NoAssocPt, meETLTrackMatchedTPmtd1Pt, meETLTPmtd1NoAssocPtMatchEff_);
  
  
  
  meETLTPmtd2CorrectAssocEtaMatchEff_ = ibook.book1D("ETLTPmtd2CorrectAssocEtaMatchEff",
						     "Track efficiency TP-mtd hit (2 sim hits), correct reco match VS Eta",
						     meETLTrackMatchedTPmtd2Eta->getNbinsX(),
						     meETLTrackMatchedTPmtd2Eta->getTH1()->GetXaxis()->GetXmin(),
						     meETLTrackMatchedTPmtd2Eta->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd2CorrectAssocEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd2CorrectAssocEta, meETLTrackMatchedTPmtd2Eta, meETLTPmtd2CorrectAssocEtaMatchEff_);
  
  meETLTPmtd2CorrectAssocPtMatchEff_ = ibook.book1D("ETLTPmtd2CorrectAssocPtMatchEff",
						    "Track efficiency TP-mtd hit (2 sim hits), correct reco match VS Pt",
						    meETLTrackMatchedTPmtd2Pt->getNbinsX(),
						    meETLTrackMatchedTPmtd2Pt->getTH1()->GetXaxis()->GetXmin(),
						    meETLTrackMatchedTPmtd2Pt->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd2CorrectAssocPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd2CorrectAssocPt, meETLTrackMatchedTPmtd2Pt, meETLTPmtd2CorrectAssocPtMatchEff_);
  
  meETLTPmtd2WrongAssocEtaMatchEff_ = ibook.book1D("ETLTPmtd2WrongAssocEtaMatchEff",
						   "Track efficiency TP-mtd hit (2 sim hits), incorrect reco match VS Eta",
						   meETLTrackMatchedTPmtd2Eta->getNbinsX(),
						   meETLTrackMatchedTPmtd2Eta->getTH1()->GetXaxis()->GetXmin(),
						   meETLTrackMatchedTPmtd2Eta->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd2WrongAssocEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd2WrongAssocEta, meETLTrackMatchedTPmtd2Eta, meETLTPmtd2WrongAssocEtaMatchEff_);
  
  meETLTPmtd2WrongAssocPtMatchEff_ = ibook.book1D("ETLTPmtd2WrongAssocPtMatchEff",
						  "Track efficiency TP-mtd hit (2 sim hits), incorrect reco match VS Pt",
						  meETLTrackMatchedTPmtd2Pt->getNbinsX(),
						  meETLTrackMatchedTPmtd2Pt->getTH1()->GetXaxis()->GetXmin(),
						  meETLTrackMatchedTPmtd2Pt->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd2WrongAssocPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd2WrongAssocPt, meETLTrackMatchedTPmtd2Pt, meETLTPmtd2WrongAssocPtMatchEff_);
  
  
  meETLTPmtd2NoAssocEtaMatchEff_ = ibook.book1D("ETLTPmtd2NoAssocEtaMatchEff",
						"Track efficiency TP-mtd hit (2 sim hits), no reco match VS Eta",
						meETLTrackMatchedTPmtd2Eta->getNbinsX(),
						meETLTrackMatchedTPmtd2Eta->getTH1()->GetXaxis()->GetXmin(),
						meETLTrackMatchedTPmtd2Eta->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd2NoAssocEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd2NoAssocEta, meETLTrackMatchedTPmtd2Eta, meETLTPmtd2NoAssocEtaMatchEff_);
  
  meETLTPmtd2NoAssocPtMatchEff_ = ibook.book1D("ETLTPmtd2NoAssocPtMatchEff",
					       "Track efficiency TP-mtd hit (2 sim hits), no reco match VS Pt",
					       meETLTrackMatchedTPmtd2Pt->getNbinsX(),
					       meETLTrackMatchedTPmtd2Pt->getTH1()->GetXaxis()->GetXmin(),
					       meETLTrackMatchedTPmtd2Pt->getTH1()->GetXaxis()->GetXmax());
  meETLTPmtd2NoAssocPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPmtd2NoAssocPt, meETLTrackMatchedTPmtd2Pt, meETLTPmtd2NoAssocPtMatchEff_);
  
  
  meETLTPnomtdEtaMatchEff_ = ibook.book1D("ETLTPnomtdEtaMatchEff",
					  "Track efficiency TP- no mtd hit, with reco match VS Eta",
					  meETLTrackMatchedTPnomtdEta->getNbinsX(),
					  meETLTrackMatchedTPnomtdEta->getTH1()->GetXaxis()->GetXmin(),
					  meETLTrackMatchedTPnomtdEta->getTH1()->GetXaxis()->GetXmax());
  meETLTPnomtdEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPnomtdAssocEta, meETLTrackMatchedTPnomtdEta, meETLTPnomtdEtaMatchEff_);
  
  meETLTPnomtdPtMatchEff_ = ibook.book1D("ETLTPnomtdPtMatchEff",
					 "Track efficiency TP- no mtd hit, with reco match VS Pt",
					 meETLTrackMatchedTPnomtdPt->getNbinsX(),
					 meETLTrackMatchedTPnomtdPt->getTH1()->GetXaxis()->GetXmin(),
					 meETLTrackMatchedTPnomtdPt->getTH1()->GetXaxis()->GetXmax());
  meETLTPnomtdPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackMatchedTPnomtdAssocPt, meETLTrackMatchedTPnomtdPt, meETLTPnomtdPtMatchEff_);
  //
    
  meNoTimeFraction_ = ibook.book1D("NoTimeFraction",
                                   "Fraction of tracks with MTD hits and no time associated; Num. of hits",
                                   meTrackNumHits->getNbinsX(),
                                   meTrackNumHits->getTH1()->GetXaxis()->GetXmin(),
                                   meTrackNumHits->getTH1()->GetXaxis()->GetXmax());
  meNoTimeFraction_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackNumHitsNT, meTrackNumHits, meNoTimeFraction_);

  meBtlEtaEff_->getTH1()->SetMinimum(0.);
  meBtlPhiEff_->getTH1()->SetMinimum(0.);
  meBtlPtEff_->getTH1()->SetMinimum(0.);
  meEtlEtaEff_->getTH1()->SetMinimum(0.);
  meEtlPhiEff_->getTH1()->SetMinimum(0.);
  meEtlPtEff_->getTH1()->SetMinimum(0.);
  meEtlEtaEff2_->getTH1()->SetMinimum(0.);
  meEtlPhiEff2_->getTH1()->SetMinimum(0.);
  meEtlPtEff2_->getTH1()->SetMinimum(0.);
  for (int i = 0; i < 2; i++) {
    meEtlEtaEffLowPt_[i]->getTH1()->SetMinimum(0.);
    meEtlEtaEff2LowPt_[i]->getTH1()->SetMinimum(0.);
  }

  meExtraPhiAtBTLEff_ = ibook.book1D("ExtraPhiAtBTLEff",
                                     "Efficiency to match hits at BTL surface of extrapolated tracks associated to LV",
                                     meExtraPhiAtBTL->getNbinsX(),
                                     meExtraPhiAtBTL->getTH1()->GetXaxis()->GetXmin(),
                                     meExtraPhiAtBTL->getTH1()->GetXaxis()->GetXmax());
  meExtraPhiAtBTLEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraPhiAtBTLmatched, meExtraPhiAtBTL, meExtraPhiAtBTLEff_);

  normalize(meExtraBTLeneInCone, 1.);

  meExtraMTDfailExtenderEtaEff_ =
      ibook.book1D("ExtraMTDfailExtenderEtaEff",
                   "Track associated to LV extrapolated at MTD surface no extender efficiency VS Eta;Eta;Efficiency",
                   meTrackMatchedTPEffEtaTotLV->getNbinsX(),
                   meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmin(),
                   meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmax());
  meExtraMTDfailExtenderEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraMTDfailExtenderEta, meTrackMatchedTPEffEtaTotLV, meExtraMTDfailExtenderEtaEff_);

  meExtraMTDfailExtenderPtEff_ = ibook.book1D(
      "ExtraMTDfailExtenderPtEff",
      "Track associated to LV extrapolated at MTD surface no extender efficiency VS Pt;Pt [GeV];Efficiency",
      meTrackMatchedTPEffPtTotLV->getNbinsX(),
      meTrackMatchedTPEffPtTotLV->getTH1()->GetXaxis()->GetXmin(),
      meTrackMatchedTPEffPtTotLV->getTH1()->GetXaxis()->GetXmax());
  meExtraMTDfailExtenderPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraMTDfailExtenderPt, meTrackMatchedTPEffPtTotLV, meExtraMTDfailExtenderPtEff_);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ----------
void MtdTracksHarvester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/Tracks/");

  descriptions.add("MtdTracksPostProcessor", desc);
}

DEFINE_FWK_MODULE(MtdTracksHarvester);
