#ifndef Bianchi_TutorialDAS2012_MuTauAnalyzer_h
#define Bianchi_TutorialDAS2012_MuTauAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/LorentzVector.h"



#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <TRandom3.h>


#include <string>
#include <utility>
#include <map>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;

class MuTauAnalyzer : public edm::EDAnalyzer{


 public:

  struct more {
    bool operator() (const double& lhs, const double& rhs) const
    {return lhs>rhs;}
  };

  explicit MuTauAnalyzer(const edm::ParameterSet&);
  ~MuTauAnalyzer();

  unsigned int jetID( const pat::Jet* jet );

  void beginJob() ;
  void analyze(const edm::Event&  iEvent, const edm::EventSetup& iSetup);
  void endJob() ;

 private:

  TFile* file_;
  TTree* tree_;

  edm::InputTag diTauTag_;
  edm::InputTag jetsTag_;
  edm::InputTag newJetsTag_;
  edm::InputTag rawMetTag_;
  edm::InputTag muonsTag_;
  edm::InputTag muonsRelTag_;
  edm::InputTag verticesTag_;
  edm::InputTag triggerResultsTag_;

  bool isMC_;
  bool verbose_;
  float minCorrPt_;
  float minJetID_;
  float deltaRLegJet_;

  // trigger-matching bits, one entry per trigger filter; the sequence of paths
  // is hardcoded in the .cc; 0=FAIL, 1=PASS
  std::vector< int >* tauXTriggers_;
  // trigger bits, one entry per HLT path; the sequence of paths
  // is hardcoded in the .cc; 0=FAIL, 1=PASS
  std::vector< int >* triggerBits_;

  // vector of four-momenta of reconstructed jets, ordered by decreasing pt
  std::vector< LV >* jetsIDP4_;
  // vector of four-momenta of generated jets matched to the reco jets
  std::vector< LV >* genJetsIDP4_; 
  // vector containing the visible four-momentum of the selected di-tau
  std::vector< LV >* diTauVisP4_; 
  // vector containing the full four-momentum of the selected di-tau
  std::vector< LV >* diTauSVfitP4_; 
  // vector containing the visible muon four-momentum (element [0]) and 
  // tau four-momentum (element [1]) 
  std::vector< LV >* diTauLegsP4_; 
  // vector containing the visible generor-level muon four-momentum (element [0]) and 
  // the generor-level tau four-momentum (element [1]) 
  std::vector< LV >* genDiTauLegsP4_;
  // vector containing the full four-momentum of the missing transverse energy (MET)
  std::vector< LV >* METP4_;
  // vector containing the sum of all neutrinos four-momentum
  std::vector< LV >* genMETP4_;
  // vector containing the four-momentum of the heavy generator-level boson (Z,W or H)
  std::vector< LV >* genVP4_;

  // an integer saying the pdgId of the heavy boson TIMES the pdgId of the lepton it decays to
  int genDecay_;
  // run, event and lumi numbers to identify an event
  unsigned long run_,event_,lumi_;
  // scalar SumPt of all PF candidates around the muon, of type 
  // PFCandidate::h, PFCandidate::nh and PFCandidate::gamma 
  float chIsoLeg1v2_,nhIsoLeg1v2_,phIsoLeg1v2_;
  // scalar SumPt of all PF candidates around the muon of type 
  // PFCandidate::h FROM PILEUP VERTICES!!! 
  float nhIsoPULeg1v2_;
  // transverse mass between the muon and the MET
  float MtLeg1_;
  // number of reconstructed primary vertices
  float numPV_;
  // number of di-taus (used for bookkeeping of the combinatorics)
  int numOfLooseIsoDiTaus_;
  // reconstructed decay mode of the tau: 1prong (0), 1prong+pi0s (1) or 3prongs (2)
  int decayMode_;
  // visible mass of the reconstructed tau lepton
  float visibleTauMass_;
  // number of charged particles inside the reconstructed tau (1 or 3)
  int signalPFChargedHadrCands_;
  // number of photons inside the reconstructed tau (0,1,2,..)
  int signalPFGammaCands_;
  // tightest of the tau isolation working points passed by the tau candidate (0=LOOSE,1=MEDIUM,2=TIGHT)
  int tightestHPSDBWP_;
  // is the reconstructed tau matched to a generator-level hadronically-decaying tau ???
  int isTauLegMatched_;
  // is the reconstructed mu matched to a generator-level muon ???
  int isMuLegMatched_;
  // flag the event according to the nnumber of loose muons in the event: 0= [only one], 1= [more than one]  
  int muFlag_;
  // muon charge + tau charge
  float diTauCharge_;
  // number of actual in-time pile-up interactions simulated on top of the hard event
  int nPUVertices_;

};


#endif
