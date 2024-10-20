#include "GeneratorInterface/GenFilters/interface/VBFGenJetFilter.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include <HepMC/GenVertex.h>

// ROOT includes
#include "TMath.h"

// C++ includes
#include <iostream>

using namespace edm;
using namespace std;


VBFGenJetFilter::VBFGenJetFilter(const edm::ParameterSet& iConfig) :
oppositeHemisphere    (iConfig.getUntrackedParameter<bool>  ("oppositeHemisphere",false)),
leadJetsNoLepMass     (iConfig.getUntrackedParameter<bool>  ("leadJetsNoLepMass", false)),
ptMin                 (iConfig.getUntrackedParameter<double>("minPt",                20)),
etaMin                (iConfig.getUntrackedParameter<double>("minEta",             -5.0)),
etaMax                (iConfig.getUntrackedParameter<double>("maxEta",              5.0)),
minInvMass            (iConfig.getUntrackedParameter<double>("minInvMass",          0.0)),
maxInvMass            (iConfig.getUntrackedParameter<double>("maxInvMass",      99999.0)),
minLeadingJetsInvMass (iConfig.getUntrackedParameter<double>("minLeadingJetsInvMass",          0.0)),
maxLeadingJetsInvMass (iConfig.getUntrackedParameter<double>("maxLeadingJetsInvMass",      99999.0)),
deltaRNoLep           (iConfig.getUntrackedParameter<double>("deltaRNoLep",        0.3)),
minDeltaPhi           (iConfig.getUntrackedParameter<double>("minDeltaPhi",        -1.0)),
maxDeltaPhi           (iConfig.getUntrackedParameter<double>("maxDeltaPhi",     99999.0)),
minDeltaEta           (iConfig.getUntrackedParameter<double>("minDeltaEta",        -1.0)),
maxDeltaEta           (iConfig.getUntrackedParameter<double>("maxDeltaEta",     99999.0))
{
  
  m_inputTag_GenJetCollection       = consumes<reco::GenJetCollection>(iConfig.getUntrackedParameter<edm::InputTag>("inputTag_GenJetCollection",edm::InputTag("ak5GenJetsNoNu")));
  if(leadJetsNoLepMass) m_inputTag_GenParticleCollection  = consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genParticles",edm::InputTag("genParticles")));
//   m_inputTag_GenParticleCollection  = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));

}

VBFGenJetFilter::~VBFGenJetFilter(){
  
}


bool VBFGenJetFilter::isHardProcess(const reco::GenParticle &p)  {
  
  //status 3 in pythia6 means hard process;
  if (p.status()==3) return true;
  
  //hard process codes for pythia8 are 21-29 inclusive (currently 21,22,23,24 are used)
  if (p.status()>20 && p.status()<30) return true;
  
  //if this is a final state or decayed particle,
  //check if direct mother is a resonance decay in pythia8 but exclude FSR branchings
  //(In pythia8 if a resonance decay product did not undergo any further branchings
  //it will be directly stored as status 1 or 2 without any status 23 copy)
  if (p.status()==1 || p.status()==2) {
//     const reco::GenParticle *um = mother(p);
    const reco::GenParticle *um = static_cast<const reco::GenParticle*>(p.mother(0));
    if (um) {
      const reco::GenParticle *firstcopy = firstCopy(*um);
      bool fromResonance = firstcopy && firstcopy->status()==22;
      
      const reco::GenParticle *umNext = nextCopy(*um);
      bool fsrBranching = umNext && umNext->status()>50 && umNext->status()<60;
      
      if (fromResonance && !fsrBranching) return true;
    }
  }
  
  return false;
  
}

/////////////////////////////////////////////////////////////////////////////

const reco::GenParticle * VBFGenJetFilter::firstCopy(const reco::GenParticle &p)  {
  const reco::GenParticle *pcopy = &p;
  std::unordered_set<const reco::GenParticle*> dupCheck;
  while (previousCopy(*pcopy)) {
    dupCheck.insert(pcopy);
    pcopy = previousCopy(*pcopy);
    if (dupCheck.count(pcopy)) return 0;
  }
  return pcopy;    
}

/////////////////////////////////////////////////////////////////////////////

const reco::GenParticle * VBFGenJetFilter::previousCopy(const reco::GenParticle &p)  {
  
  const unsigned int nmoth = p.numberOfMothers();
  for (unsigned int imoth = 0; imoth<nmoth; ++imoth) {
    const reco::GenParticle *moth = static_cast<const reco::GenParticle*>(p.mother(imoth));//mother(p,imoth);
    if (moth->pdgId()==p.pdgId()) {
      return moth;
    }
  }
  
  return 0;     
}   

/////////////////////////////////////////////////////////////////////////////

const reco::GenParticle * VBFGenJetFilter::nextCopy(const reco::GenParticle &p)  {
  
  const unsigned int ndau = p.numberOfDaughters();
  for (unsigned int idau = 0; idau<ndau; ++idau) {
    const reco::GenParticle *dau = static_cast<const reco::GenParticle*>(p.daughter(idau));//daughter(p,idau);
    if (dau->pdgId()==p.pdgId()) {
      return dau;
    }
  }
  
  return 0;     
} 

/////////////////////////////////////////////////////////////////////////////

vector<const reco::GenParticle*> VBFGenJetFilter::filterGenLeptons(const vector<reco::GenParticle>* particles){
  vector<const reco::GenParticle*> out;
  
  



  for(const auto & p : *particles){
      
      int absPdgId = std::abs(p.pdgId());
      
      if(((absPdgId == 11) || (absPdgId == 13) || (absPdgId == 15)) && isHardProcess(p)) {
          out.push_back(&p);              
      }
       
          
  }     
  return out;
}


/////////////////////////////////////////////////////////////////////////////


vector<const reco::GenJet*> VBFGenJetFilter::filterGenJets(const vector<reco::GenJet>* jets){
  
  vector<const reco::GenJet*> out;
  
  for(unsigned i=0; i<jets->size(); i++){
    
    const reco::GenJet* j = &((*jets)[i]);
    
    if(j->p4().pt() >ptMin &&  j->p4().eta()>etaMin && j->p4().eta()<etaMax)
    {
      out.push_back(j);
    }
  }
  
  return out;
}


/////////////////////////////////////////////////////////////////////////////




// ------------ method called to skim the data  ------------
bool VBFGenJetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  Handle< vector<reco::GenJet> > handleGenJets;
  iEvent.getByToken(m_inputTag_GenJetCollection, handleGenJets);
  const vector<reco::GenJet>* genJets = handleGenJets.product();
  
        
        
  // Getting filtered generator jets
  vector<const reco::GenJet*> filGenJets = filterGenJets(genJets);
  
  // If we do not find at least 2 jets veto the event
  if(filGenJets.size()<2){return false;}
  
  
  // Testing dijet mass   
  if(leadJetsNoLepMass) { 
        
    Handle<reco::GenParticleCollection> genParticelesCollection;
    iEvent.getByToken(m_inputTag_GenParticleCollection, genParticelesCollection);
    const vector<reco::GenParticle>*  genParticles = genParticelesCollection.product();
            
    
      // Getting filtered generator muons
      vector<const reco::GenParticle*> filGenLep = filterGenLeptons(genParticles);

      
      // Getting p4 of jet with no lepton
      vector<math::XYZTLorentzVector> genJetsWithoutLeptonsP4;
      unsigned int jetIdx = 0;
      
      
      

      while(genJetsWithoutLeptonsP4.size()<2 && jetIdx < filGenJets.size()) {
          bool jetWhitoutLep = true;
          const math::XYZTLorentzVector & p4J= (filGenJets[jetIdx])->p4();
          for(unsigned int i = 0; i < filGenLep.size() && jetWhitoutLep; ++i) {
              if(reco::deltaR2((filGenLep[i])->p4(), p4J) < deltaRNoLep*deltaRNoLep)
              
                  jetWhitoutLep = false;
          }
          
          if (jetWhitoutLep)  genJetsWithoutLeptonsP4.push_back(p4J);
          ++jetIdx;
      }
      
      // Checking the invariant mass of the leading jets
      if (genJetsWithoutLeptonsP4.size() < 2) return false;
      float invMassLeadingJet = (genJetsWithoutLeptonsP4[0] + genJetsWithoutLeptonsP4[1]).M();
      if ( invMassLeadingJet > minLeadingJetsInvMass  && invMassLeadingJet < maxLeadingJetsInvMass) return true;
      else return false;
      
  }

  

  
  for(unsigned a=0; a<filGenJets.size(); a++){
    for(unsigned b=a+1; b<filGenJets.size(); b++){    
      
      const reco::GenJet* pA = filGenJets[a];
      const reco::GenJet* pB = filGenJets[b];
      
      // Getting the dijet vector
      math::XYZTLorentzVector diJet = pA->p4() + pB->p4();
      
      // Testing opposite hemispheres
      double dijetProd = pA->p4().eta()*pB->p4().eta();
      if(oppositeHemisphere && dijetProd>=0){continue;}
      
      // Testing dijet mass
      double invMass = diJet.mass();
      if(invMass<=minInvMass || invMass>maxInvMass){continue;}
      
      
      // Testing dijet delta eta
      double dEta = fabs(pA->p4().eta()-pB->p4().eta());
      if(dEta<=minDeltaEta || dEta>maxDeltaEta){continue;}

      // Testing dijet delta phi
      double dPhi = fabs(reco::deltaPhi(pA->p4().phi(),pB->p4().phi()));
      if(dPhi<=minDeltaPhi || dPhi>maxDeltaPhi){continue;}
      
      return true;
    }
  }

  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(VBFGenJetFilter);
