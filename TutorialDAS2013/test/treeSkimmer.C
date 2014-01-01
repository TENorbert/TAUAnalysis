#include "TTree.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TBenchmark.h"
#include "TGraph.h"
#include "TVectorT.h"
#include "TMultiGraph.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TLorentzVector.h"

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include <vector>
#include <utility>
#include <map>
#include <algorithm>

#define SAVE   true
#define MINPt1 20.0 
#define MINPt2 20.0
#define PtVETO 30.0
#define MAXEta  4.5 

using namespace ROOT::Math;
using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;


float pileupWeight( int intimepileup_ ){

  
 
  std::vector< float > Data2011LumiExt;
  Double_t Data2011LumiExt_f[50] = {
    0.00563937,
    0.0229973, 
    0.0517398, 
    0.0847196, 
    0.112657, 
    0.128767, 
    0.130759, 
    0.120404, 
    0.101874, 
    0.0799242,
    0.0585292,
    0.0402185,
    0.0260466,
    0.0159597,
    0.00928456,
    0.00514434,
    0.00272259,
    0.00137995,
    0.000671453, 
    0.000314329, 
    0.000141851, 
    6.18205e-05, 
    2.6061e-05, 
    1.06425e-05,
    4.21557e-06,
    1.62162e-06,
    6.06438e-07,
    2.2069e-07, 
    7.82201e-08,
    2.7023e-08, 
    9.10621e-09,
    2.99511e-09,
    9.62094e-10,
    3.01987e-10,
    1.31916e-10,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };

  std::vector< float > Fall11Lumi ;
  Double_t Fall11Lumi_f[50] = {
    0.003388501,
    0.010357558,
    0.024724258,
    0.042348605,
    0.058279812,
    0.068851751,
    0.072914824,
    0.071579609,
    0.066811668,
    0.060672356,
    0.054528356,
    0.04919354,
    0.044886042,
    0.041341896,
    0.0384679,
    0.035871463,
    0.03341952,
    0.030915649,
    0.028395374,
    0.025798107,
    0.023237445,
    0.020602754,
    0.0180688,
    0.015559693,
    0.013211063,
    0.010964293,
    0.008920993,
    0.007080504,
    0.005499239,
    0.004187022,
    0.003096474,
    0.002237361,
    0.001566428,
    0.001074149,
    0.000721755,
    0.000470838,
    0.00030268,
    0.000184665,
    0.000112883,
    6.74043E-05,
    3.82178E-05,
    2.22847E-05,
    1.20933E-05,
    6.96173E-06,
    3.4689E-06,
    1.96172E-06,
    8.49283E-07,
    5.02393E-07,
    2.15311E-07,
    9.56938E-08
  };

 
  for( int i=0; i<50; i++) {
    Data2011LumiExt.push_back(Data2011LumiExt_f[i]);
    Fall11Lumi.push_back(Fall11Lumi_f[i]);
  }

  int intimepileup = intimepileup_>49 ? 49 : intimepileup_;

  return Data2011LumiExt[intimepileup]/Fall11Lumi[intimepileup];

}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////



void makeTrees_MuTau( int index = 0 ){
 
  cout <<  "Now analyzing sample with index " << index << endl;

  TFile corrections("../data/corrections.root");
  TF1 *ratioMuIDIsoBL     = (TF1*)corrections.Get("ratioMuIDIsoBL");
  TF1 *ratioMuIDIsoEC     = (TF1*)corrections.Get("ratioMuIDIsoEC");
  TF1 *ratioMuHLTBL       = (TF1*)corrections.Get("ratioMuHLTBL");
  TF1 *ratioMuHLTEC       = (TF1*)corrections.Get("ratioMuHLTEC");
  TF1 *ratioTauHLTBL      = (TF1*)corrections.Get("ratioTauHLTBL");  
  TF1 *ratioTauHLTEC      = (TF1*)corrections.Get("ratioTauHLTEC");  

  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;

  std::vector<std::string> samples;
  std::vector<float> crossSec;

  // samples & x-sections & skim1 
  samples.push_back("Data");        crossSec.push_back(     0  );                          
  samples.push_back("DYJets");      crossSec.push_back(  3048.0     * 0.00963); 
  samples.push_back("TTJets");      crossSec.push_back(   157.5     * 0.0210);  
  samples.push_back("WJets");       crossSec.push_back( 31314.0     * 0.00126);   
  samples.push_back("VBFH130");     crossSec.push_back(     0.0632  * 0.066);   
  samples.push_back("GGFH130");     crossSec.push_back(     0.774   * 0.055);  
  samples.push_back("QCD");         crossSec.push_back(   641              );
 
  // normalization Lumi
  Float_t Lumi= 100;

  string currentSample = samples[index];
  TString sample(currentSample.c_str());

  TString outName = "treeSkimmedMuTau_"+sample+".root";
  TFile *outFile = new TFile(outName,"UPDATE");

  TTree* outTreePtOrd     = new TTree("outTreePtOrd","tree");

  // kinematical variables of leading and trailing jet (ordered by pt)
  // pt, eta, |eta1-eta2|, mass of the di-jet system, azimuthal angle, phi1-phi2  
  float pt1,pt2,eta1,eta2,Deta,Mjj,Dphi,phi1,phi2;
  // number of jets above 30 or 20 GeV
  int nJets30, nJets20;

  // pt, eta and phi of the third veto
  float ptVeto, etaVeto, phiVeto; 
  // bool saying if the third jet lies within the leading and trailing jets
  int isVetoInJets;

  // diTau related variables
  // full mass of the di-tau
  float diTauSVFitMass;
  // visible mass, pt , eta and phi of the di-tau
  float diTauVisMass,diTauVisPt,diTauVisEta,diTauVisPhi;

  // taus/MET related variables
  // pt and eta of the muon (1) and tau (2)
  float ptL1,ptL2,etaL1,etaL2,phiL1,phiL2,dPhiL1L2;
  //  As in the input TTree..
  float diTauCharge_;
  //  As in the input TTree..
  float MtLeg1_;
  // MET pt
  float MEt;
  // relative isolation variable: w/o and w/ pile-up correction
  float combRelIsoLeg1DBeta;
  //  As in the input TTree..
  int tightestHPSDBWP_, decayMode_;
  float visibleTauMass_;

  // event-related variables
  // numnber of reconstructed vertices
  float numPV_ ;
  // weight to normalize the sample to 100pb-1
  float sampleWeight;
  // pile-up reweighting 
  float puWeight;
  // Higgs boson pt spectrum reweighting
  float HqTWeight;

  // object-related weights and triggers
  // trigger bit: 1=FIRED ; 0=NOT FIRED
  float HLTx;
  // trigger matching: 1=MATCHED; 0=UNMATCHED
  float HLTmatch;
  // Trigger turn-on correction
  float HLTweightTau, HLTweightMu;
  // Muon and tau isolation(X)identification data/MC scale factor
  float SFMu, SFTau;
  // As in the input TTree..
  int muFlag_;
  // As in the input TTree..
  int genDecay_;

  // As in the input TTree..
  ULong64_t event_,run_,lumi_;

  outTreePtOrd->Branch("pt1",  &pt1,"pt1/F");
  outTreePtOrd->Branch("pt2",  &pt2,"pt2/F");
  outTreePtOrd->Branch("eta1", &eta1,"eta1/F");
  outTreePtOrd->Branch("eta2", &eta2,"eta2/F");
  outTreePtOrd->Branch("phi1", &phi1,"phi1/F");
  outTreePtOrd->Branch("phi2", &phi2,"phi2/F");
  outTreePtOrd->Branch("Deta", &Deta,"Deta/F");
  outTreePtOrd->Branch("Dphi", &Dphi,"Dphi/F");
  outTreePtOrd->Branch("Mjj",  &Mjj,"Mjj/F");

  outTreePtOrd->Branch("nJets30",       &nJets30,  "nJets30/I");  
  outTreePtOrd->Branch("nJets20",       &nJets20,  "nJets20/I");  

  outTreePtOrd->Branch("ptVeto",  &ptVeto, "ptVeto/F");
  outTreePtOrd->Branch("phiVeto", &phiVeto,"phiVeto/F");
  outTreePtOrd->Branch("etaVeto", &etaVeto,"etaVeto/F");
  outTreePtOrd->Branch("isVetoInJets", &isVetoInJets,"isVetoInJets/I");
 
  outTreePtOrd->Branch("diTauSVFitMass",&diTauSVFitMass,"diTauSVFitMass/F");
  
  outTreePtOrd->Branch("diTauVisMass",&diTauVisMass,"diTauVisMass/F");
  outTreePtOrd->Branch("diTauVisPt",  &diTauVisPt,"diTauVisPt/F");
  outTreePtOrd->Branch("diTauVisEta", &diTauVisEta,"diTauVisEta/F");
  outTreePtOrd->Branch("diTauVisPhi", &diTauVisPhi,"diTauVisPhi/F");

  outTreePtOrd->Branch("etaL1",   &etaL1,"etaL1/F");
  outTreePtOrd->Branch("etaL2",   &etaL2,"etaL2/F");
  outTreePtOrd->Branch("ptL1",    &ptL1,"ptL1/F");
  outTreePtOrd->Branch("ptL2",    &ptL2,"ptL2/F");
  outTreePtOrd->Branch("phiL1",   &phiL1,"phiL1/F");
  outTreePtOrd->Branch("phiL2",   &phiL2,"phiL2/F");
  outTreePtOrd->Branch("dPhiL1L2",&dPhiL1L2,"dPhiL1L2/F");

  outTreePtOrd->Branch("diTauCharge", &diTauCharge_,"diTauCharge/F");
  outTreePtOrd->Branch("MtLeg1",      &MtLeg1_,"MtLeg1/F");
  outTreePtOrd->Branch("MEt",         &MEt,"MEt/F");

  outTreePtOrd->Branch("combRelIsoLeg1DBeta",&combRelIsoLeg1DBeta,"combRelIsoLeg1DBeta/F");

  outTreePtOrd->Branch("tightestHPSDBWP",  &tightestHPSDBWP_,"tightestHPSDBWP/I");
  outTreePtOrd->Branch("decayMode",        &decayMode_,"decayMode/I");
  outTreePtOrd->Branch("visibleTauMass",&visibleTauMass_,"visibleTauMass/F");

  outTreePtOrd->Branch("numPV",              &numPV_,"numPV/F");
  outTreePtOrd->Branch("sampleWeight",       &sampleWeight,"sampleWeight/F"); 
  outTreePtOrd->Branch("puWeight",           &puWeight,"puWeight/F");
  outTreePtOrd->Branch("HqTWeight",          &HqTWeight,"HqTWeight/F");

  outTreePtOrd->Branch("HLTx",         &HLTx, "HLTx/F");
  outTreePtOrd->Branch("HLTmatch",     &HLTmatch,"HLTmatch/F");
  outTreePtOrd->Branch("HLTweightMu",  &HLTweightMu,"HLTweightMu/F");
  outTreePtOrd->Branch("HLTweightTau", &HLTweightTau,"HLTweightTau/F");
  outTreePtOrd->Branch("SFTau",        &SFTau,"SFTau/F");
  outTreePtOrd->Branch("SFMu",         &SFMu,"SFMu/F");

  outTreePtOrd->Branch("muFlag",          &muFlag_,"muFlag/I"); 
  outTreePtOrd->Branch("genDecay",        &genDecay_,"genDecay/I");

  outTreePtOrd->Branch("event",&event_,"event/l");
  outTreePtOrd->Branch("run",  &run_,  "run/l");
  outTreePtOrd->Branch("lumi", &lumi_, "lumi/l");
 
  string currentInName = "./treeMuTau_"+samples[index]+".root" ;

  TString inName(currentInName.c_str());
  TFile* file   = new TFile(inName,"READ");
  if(file->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }

  TString treeName("muTauStreamAnalyzerRaw/tree");

  TTree* currentTree = (TTree*)file->Get(treeName);
  int nEntries    = currentTree->GetEntries() ;
  int nEventsRead = ((TH1F*)file->Get("allEventsFilter/totalEvents"))->GetBinContent(1) ;
  float crossSection = crossSec[index] ; 
  float scaleFactor = (crossSection != 0) ? Lumi / (  float(nEventsRead)/crossSection )  : 1.0;

  cout << "Processing sample " << sample.Data() << endl;
  cout<< "nEventsRead = " << nEventsRead << endl;
  cout<< "nEntries    = " << nEntries << endl;
  cout<< "crossSection " << crossSection << " pb ==> scaleFactor " << scaleFactor << endl;

  // jets
  currentTree->SetBranchStatus("jetsIDP4"              ,1);
  currentTree->SetBranchStatus("genJetsIDP4"           ,0);

  // diTaus
  currentTree->SetBranchStatus("diTauVisP4"            ,1);
  currentTree->SetBranchStatus("diTauSVfitP4"          ,1);
  currentTree->SetBranchStatus("diTauCharge"           ,1);

  // taus
  currentTree->SetBranchStatus("diTauLegsP4"           ,1);
  currentTree->SetBranchStatus("genDiTauLegsP4"        ,1);
  currentTree->SetBranchStatus("chIsoLeg1v2"           ,1);
  currentTree->SetBranchStatus("nhIsoLeg1v2"           ,1);
  currentTree->SetBranchStatus("phIsoLeg1v2"           ,1);
  currentTree->SetBranchStatus("nhIsoPULeg1v2"         ,1);
  currentTree->SetBranchStatus("decayMode"             ,1);
  currentTree->SetBranchStatus("tightestHPSDBWP"       ,1);
  currentTree->SetBranchStatus("visibleTauMass"        ,1);
  currentTree->SetBranchStatus("isTauLegMatched"       ,0);
  currentTree->SetBranchStatus("isMuLegMatched"        ,0);
  currentTree->SetBranchStatus("signalPFChargedHadrCands"    ,1);
  currentTree->SetBranchStatus("signalPFGammaCands"    ,1);

  // MET
  currentTree->SetBranchStatus("METP4"                 ,1);
  currentTree->SetBranchStatus("genMETP4"              ,1);
  currentTree->SetBranchStatus("MtLeg1"                ,1);

  // generator-level boson
  currentTree->SetBranchStatus("genVP4"                ,1);
  currentTree->SetBranchStatus("genDecay"              ,1);

  // extra-muons
  currentTree->SetBranchStatus("muFlag"                ,1);

  // event-dependent variables
  currentTree->SetBranchStatus("numPV"                 ,1);
  currentTree->SetBranchStatus("nPUVertices"           ,1);
  currentTree->SetBranchStatus("numOfLooseIsoDiTaus"   ,0);
  currentTree->SetBranchStatus("run"                   ,1);
  currentTree->SetBranchStatus("event"                 ,1);
  currentTree->SetBranchStatus("lumi"                  ,1);

  // triggers
  currentTree->SetBranchStatus("tauXTriggers"          ,1);
  currentTree->SetBranchStatus("triggerBits"           ,1);

  ////////////////////////////////////////////////////////////////////

  std::vector< LV >* jets           = new std::vector< LV >();
  currentTree->SetBranchAddress("jetsIDP4",          &jets);

  std::vector< LV >* diTauLegsP4    = new std::vector< LV >();
  currentTree->SetBranchAddress("diTauLegsP4",     &diTauLegsP4);

  std::vector< LV >* diTauVisP4     = new std::vector< LV >();
  currentTree->SetBranchAddress("diTauVisP4",      &diTauVisP4);

  std::vector< LV >* diTauSVfitP4     = new std::vector< LV >();
  currentTree->SetBranchAddress("diTauSVfitP4",    &diTauSVfitP4);

  std::vector< LV >* genDiTauLegsP4 = new std::vector< LV >();
  currentTree->SetBranchAddress("genDiTauLegsP4",    &genDiTauLegsP4);

  std::vector< LV >* genVP4         = new std::vector< LV >();
  currentTree->SetBranchAddress("genVP4",          &genVP4);

  std::vector< LV >* METP4          = new std::vector< LV >();
  currentTree->SetBranchAddress("METP4",           &METP4);

  std::vector< LV >* genMETP4          = new std::vector< LV >();
  currentTree->SetBranchAddress("genMETP4",        &genMETP4);

  std::vector< int >* tauXTriggers  = new std::vector< int >();
  currentTree->SetBranchAddress("tauXTriggers",    &tauXTriggers);
  std::vector< int >* triggerBits   = new std::vector< int >();
  currentTree->SetBranchAddress("triggerBits",     &triggerBits);

  // auxiliary float to store branch values
  float diTauCharge;
  int tightestHPSDBWP, decayMode;
  float numPV;
  int signalPFChargedHadrCands, signalPFGammaCands;
  int muFlag,genDecay, nPUVertices;
  float visibleTauMass;
  float chIsoLeg1,nhIsoLeg1,phIsoLeg1; 
  float nhIsoPULeg1; 
  ULong64_t event,run,lumi;

  currentTree->SetBranchAddress("chIsoLeg1v2",          &chIsoLeg1);
  currentTree->SetBranchAddress("nhIsoLeg1v2",          &nhIsoLeg1);
  currentTree->SetBranchAddress("phIsoLeg1v2",          &phIsoLeg1);
  currentTree->SetBranchAddress("nhIsoPULeg1v2",        &nhIsoPULeg1);
  currentTree->SetBranchAddress("tightestHPSDBWP",      &tightestHPSDBWP);
  currentTree->SetBranchAddress("diTauCharge",          &diTauCharge);
  currentTree->SetBranchAddress("numPV",                &numPV);
  currentTree->SetBranchAddress("event",                &event);
  currentTree->SetBranchAddress("run",                  &run);
  currentTree->SetBranchAddress("lumi",                 &lumi);
  currentTree->SetBranchAddress("nPUVertices",          &nPUVertices);
  currentTree->SetBranchAddress("genDecay",             &genDecay);
  currentTree->SetBranchAddress("decayMode",            &decayMode);
  currentTree->SetBranchAddress("muFlag",               &muFlag);
  currentTree->SetBranchAddress("visibleTauMass",       &visibleTauMass);
  currentTree->SetBranchAddress("signalPFChargedHadrCands"    ,&signalPFChargedHadrCands);
  currentTree->SetBranchAddress("signalPFGammaCands"   ,&signalPFGammaCands);

  TFile* HqT      = 0;
  int mH          = 130;
  TH1F* histo     = 0;
  if(samples[index].find("GGFH130")!=string::npos){
    mH = 130;
    cout << "Reweighting powheg with HqT mH=" << mH << endl;
    HqT = new TFile(Form("../data/weight_ptH_%d.root", mH));
    if(!HqT) cout << "Cannot find HqT file..." << endl;
    else{
      histo = (TH1F*)(HqT->Get(Form("powheg_weight/weight_hqt_fehipro_fit_%d",mH)));
    }
  }


  for (int n = 0; n < nEntries ; n++) {
    
    currentTree->GetEntry(n);
    if(n%5000==0) cout << n << endl;
    
    // initialize variables filled only in the two jet case
    pt1=-99;pt2=-99;eta1=-99,eta2=-99;Deta=-99;Dphi=-99;Mjj=-99;phi1=-99;phi2=-99;
    ptVeto = -99; phiVeto= -99; etaVeto= -99;isVetoInJets=-99;

    // define the relevant jet collection 
    nJets30        = 0;
    int lead  = -99;
    int trail = -99;
    int veto  = -99;
    vector<int> indexes;
    for(int l = 0 ; l < jets->size() ; l++){
      if((*jets)[l].Pt()>MINPt1 && TMath::Abs((*jets)[l].Eta())<MAXEta)
	indexes.push_back(l);
    }
    nJets20 = indexes.size();
    if(indexes.size()>0) lead  = indexes[0];  
    if(indexes.size()>1) trail = indexes[1];  
    if(indexes.size()>2) veto  = indexes[2];  

    for(int v = 0 ; v < indexes.size() ; v++){
      if( (*jets)[indexes[v]].Pt() > 30 ) nJets30++;
    }

    // first jet
    if(lead>=0){

      pt1  = (*jets)[lead].Pt();
      eta1 = (*jets)[lead].Eta();
      phi1 = (*jets)[lead].Phi();

      // second jet
      if(trail>=0){
	pt2  = (*jets)[trail].Pt();
	eta2 = (*jets)[trail].Eta();
	phi2 = (*jets)[trail].Phi();
	Deta = abs(eta1-eta2);
	Dphi = abs((*jets)[lead].Phi()-(*jets)[trail].Phi()) > TMath::Pi() ? 
	  -abs( (*jets)[lead].Phi()-(*jets)[trail].Phi() ) + 2*TMath::Pi()  :
	  abs( (*jets)[lead].Phi()-(*jets)[trail].Phi() ) ;
	Mjj  = ((*jets)[lead]+(*jets)[trail]).M();
      }
    }

    

    ptVeto  = (veto>=0) ? (*jets)[veto].Pt() : -99;
    etaVeto = (veto>=0) ? (*jets)[veto].Eta(): -99;
    phiVeto = (veto>=0) ? (*jets)[veto].Phi(): -99;
 
    isVetoInJets = 0;
    for(int l = 0 ; l < indexes.size() ; l++){
      if(lead>=0 && trail>=0 && (l!= lead && l!= trail) &&
	 (*jets)[indexes[l]].Pt()>PtVETO && ((*jets)[indexes[l]].Eta() - eta1)*((*jets)[indexes[l]].Eta() - eta2)<=0 )
	isVetoInJets = 1;
    }

    diTauSVFitMass = (*diTauSVfitP4)[0].M();
    diTauVisMass   = (*diTauVisP4)[0].M();
    diTauVisPt     = (*diTauVisP4)[0].Pt();
    diTauVisEta    = (*diTauVisP4)[0].Eta();
    diTauVisPhi    = (*diTauVisP4)[0].Phi();
    
    ptL1     = (*diTauLegsP4)[0].Pt();
    ptL2     = (*diTauLegsP4)[1].Pt();
    etaL1    = (*diTauLegsP4)[0].Eta();
    etaL2    = (*diTauLegsP4)[1].Eta();
    phiL1    = (*diTauLegsP4)[0].Phi();
    phiL2    = (*diTauLegsP4)[1].Phi();
    dPhiL1L2 =  abs((*diTauLegsP4)[0].Phi()-(*diTauLegsP4)[1].Phi()) > TMath::Pi() ? 
      -abs( (*diTauLegsP4)[0].Phi()-(*diTauLegsP4)[1].Phi() ) + 2*TMath::Pi()  :
      abs( (*diTauLegsP4)[0].Phi()-(*diTauLegsP4)[1].Phi() ) ;

    diTauCharge_    = diTauCharge;
    
    

    ////////////////////////////////////////////////////////////////////

    float scalarSumPt     = ( *diTauLegsP4)[0].Pt() + (*METP4)[0].Pt();
    float vectorSumPt     = ((*diTauLegsP4)[0] +      (*METP4)[0]).Pt() ;

    MtLeg1_     = TMath::Sqrt( scalarSumPt*scalarSumPt - vectorSumPt*vectorSumPt ) ;
    MEt     = (*METP4)[0].Et();
    combRelIsoLeg1DBeta    = (chIsoLeg1+ std::max( nhIsoLeg1+phIsoLeg1-0.5*(nhIsoPULeg1),double(0.0)))/(*diTauLegsP4)[0].Pt();

    tightestHPSDBWP_ = tightestHPSDBWP;
    visibleTauMass_  = visibleTauMass;
    decayMode_       = decayMode;
    numPV_           = numPV;
    sampleWeight     = scaleFactor; 
    puWeight         = (std::string(sample.Data())).find("Data")!=string::npos ? 1.0 : pileupWeight(nPUVertices);   

    HqTWeight = histo!=0 ? histo->GetBinContent( histo->FindBin( (*genVP4)[0].Pt() ) ) : 1.0;

    // move HLT check here, and check for any bit
    HLTx = 0;
    for( unsigned int iBit = 0; iBit < triggerBits->size(); iBit++) {
      HLTx =  float((*triggerBits)[iBit] || HLTx);
    }
    if((std::string(sample.Data())).find("Data")!=string::npos){
      
      //HLTx =  float((*triggerBits)[0]);  //HLT_IsoMu15_LooseIsoPFTau15_v8
      bool isTriggMatched = (*tauXTriggers)[0] && (*tauXTriggers)[1] ; //hltSingleMuIsoL3IsoFiltered15 && hltOverlapFilterIsoMu15IsoPFTau15
      HLTmatch = isTriggMatched ? 1.0 : 0.0;
      HLTweightTau = 1.0;
      HLTweightMu  = 1.0;
      SFTau        = 1.0;
      SFMu         = 1.0;
      
    } 
    else{

      //HLTx  =  float((*triggerBits)[0]); //HLT_IsoMu15_LooseIsoPFTau15_v9
      bool isTriggMatched = (*tauXTriggers)[0] && (*tauXTriggers)[1] ; //hltSingleMuIsoL3IsoFiltered15 && hltOverlapFilterIsoMu15IsoPFTau15
      HLTmatch = isTriggMatched ? 1.0 : 0.0;

      HLTweightTau  = TMath::Abs((*diTauLegsP4)[1].Eta())<1.5 ?  
	ratioTauHLTBL->Eval( (*diTauLegsP4)[1].Pt() ) : ratioTauHLTEC->Eval( (*diTauLegsP4)[1].Pt() );
      HLTweightMu  = TMath::Abs((*diTauLegsP4)[0].Eta())<1.5 ?  
	ratioMuHLTBL->Eval( (*diTauLegsP4)[0].Pt() ) : ratioMuHLTEC->Eval( (*diTauLegsP4)[0].Pt() );
      SFTau  = 1.0;
      SFMu   =  TMath::Abs((*diTauLegsP4)[0].Eta())<1.5 ?
	ratioMuIDIsoBL->Eval( (*diTauLegsP4)[0].Pt() ): 
	ratioMuIDIsoEC->Eval( (*diTauLegsP4)[0].Pt() );
      
    }
   
    muFlag_          = muFlag;
    genDecay_        = genDecay ;
    event_           = event;
    run_             = run;
    lumi_            = lumi;
    
    outTreePtOrd->Fill();
  }


  file->Close();
  
  if(SAVE) outFile->Write();
  outFile->Close();

  delete jets; delete diTauLegsP4; delete diTauVisP4; delete diTauSVfitP4; delete genDiTauLegsP4;
  delete tauXTriggers; delete triggerBits;
  delete METP4; delete genVP4; delete genMETP4;
  delete HqT;
  
  return;

}


void doAllSamples(){
 
  for( unsigned int k = 0; k < 6 ; k++) makeTrees_MuTau(k);

  return;

}

