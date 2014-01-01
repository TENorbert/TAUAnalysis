
{
  // centrally produced trees
  TString path("/uscms_data/d1/jkolb//FullSkimmedTrees/"); 
  TString treeName("outTreePtOrdRaw");
  
  // locally produced trees
  //TString path("./"); 
  //TString treeName("outTreePtOrd");

  bool verbose = false;

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptStat(0000000);
  gStyle->SetOptFit(0111);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetPalette(1);


  TFile fData(path + TString("treeSkimmedMuTau_Data.root"));
  TFile fDYJets(path + TString("treeSkimmedMuTau_DYJets.root"));
  TFile fWJets(path + TString("treeSkimmedMuTau_WJets.root"));
  TFile fTTJets(path + TString("treeSkimmedMuTau_TTJets.root"));
  TFile fVBFH130(path + TString("treeSkimmedMuTau_VBFH130.root"));
  TFile fGGFH130(path + TString("treeSkimmedMuTau_GGFH130.root"));
  if(verbose) std::cout << "Got files" << std::endl;
  TTree* tData   = (TTree*)fData.Get(treeName);
  TTree* tDYJets = (TTree*)fDYJets.Get(treeName);
  TTree* tWJets  = (TTree*)fWJets.Get(treeName);
  TTree* tTTJets = (TTree*)fTTJets.Get(treeName);
  TTree* tVBFH130= (TTree*)fVBFH130.Get(treeName);
  TTree* tGGFH130= (TTree*)fGGFH130.Get(treeName);
  if(verbose) std::cout << "Got trees" << std::endl;
  int nBins = 11;
  TArrayF bins(nBins+1);
  bins[0] = 0;    bins[1] = 30;    bins[2] = 40;    bins[3] = 50;    bins[4]  = 60;    bins[5] = 70;
  bins[6] = 80;   bins[7] = 90;    bins[8] = 100;   bins[9] = 120;   bins[10] = 150;  bins[11] = 200;
  //for(int k = 0 ; k <= nBins ; k++) bins[k] = 15*k; 
  TString variable("diTauVisMass");
  TString labels(" ; mass (GeV) ; Events");

  TH1F* hData    = new TH1F("hData"   ,labels, nBins, bins.GetArray());
  TH1F* hQCD     = new TH1F("hQCD"    ,labels, nBins, bins.GetArray());
  TH1F* hDYJets  = new TH1F("hDYJets" ,labels, nBins, bins.GetArray());
  TH1F* hWJets   = new TH1F("hWJets"  ,labels, nBins, bins.GetArray());
  TH1F* hWJetsSS = new TH1F("hWJetsSS",labels, nBins, bins.GetArray());
  TH1F* hTTJets  = new TH1F("hTTJets" ,labels, nBins, bins.GetArray());
  TH1F* hVBFH130 = new TH1F("hVBFH130",labels, nBins, bins.GetArray());
  TH1F* hGGFH130 = new TH1F("hGGFH130",labels, nBins, bins.GetArray());

  hData->SetMarkerStyle(kFullCircle);
  hQCD->SetFillColor(kMagenta-9);
  hDYJets->SetFillColor(kYellow-9);
  hWJets->SetFillColor(kRed-3);
  hTTJets->SetFillColor(kBlue-2);
  hVBFH130->SetLineWidth(3);
  hVBFH130->SetLineStyle(kDashed);
  hVBFH130->SetLineColor(kBlue);
  hGGFH130->SetLineWidth(3);
  hGGFH130->SetLineStyle(kDashed);
  hGGFH130->SetLineColor(kBlue);

  // full signal selection
  TCut sbin(         "ptL1>18 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge==0 && MtLeg1<40 && muFlag==0 && HLTx && combRelIsoLeg1DBeta<0.10)";
  // selection of opposite-sign W-enriched control region
  TCut sbinEnrichW(    "ptL1>18 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge==0 && MtLeg1>60 && muFlag==0 && HLTx && combRelIsoLeg1DBeta<0.10");
  // selection of smae-sign enriched control-region 
  TCut sbinSSEnrichW(  "ptL1>18 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge!=0 && MtLeg1>60 && muFlag==0 && HLTx && combRelIsoLeg1DBeta<0.10");
  // selection of same-sign control region
  TCut sbinSS(       "ptL1>18 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge!=0 && MtLeg1<40 && muFlag==0 && HLTx && combRelIsoLeg1DBeta<0.10");
  // selection of same sign loosely-isolated control region
  TCut sbinSSRelIso( "ptL1>18 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge!=0 && MtLeg1<40 && muFlag==0 && HLTx && combRelIsoLeg1DBeta<0.30");
 
  // renormalize the samples to the full lumi; samples are originally scaled to 1/fb
  float lumiFact    = 4700./1000;
  // opposite-sign to same-sign ratio for QCD
  float OStoSSRatio = 1.11;

  /////////////////////////////////////////////////////////////////////
  // estimation of W+jets
  std::cout << "Doing W background" << std::endl;
  TH1F* h1 = new TH1F("h1","",1,-10,10);
  tWJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbinEnrichW);
  float WsbinEnrichW  = h1->Integral()*lumiFact;
  h1->Reset();
  tWJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  float Wsbin       = h1->Integral()*lumiFact;
  h1->Reset();
  tTTJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbinEnrichW);
  float TTsbinEnrichW = h1->Integral()*lumiFact;
  h1->Reset();
  tData->Draw("etaL1>>h1",   sbinEnrichW);
  float DatasbinEnrichW = h1->Integral();
  h1->Reset();

  std::cout << "W events in sbin (MC) = " << Wsbin << endl;
  float WsbinCorr = (DatasbinEnrichW - TTsbinEnrichW)*(Wsbin/WsbinEnrichW);
  std::cout << "W events in sbin (corrected) = (" << DatasbinEnrichW << " - " << TTsbinEnrichW << " )*" << Wsbin/WsbinEnrichW << " = " << WsbinCorr << endl;
  /////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////
  // estimation of QCD
  std::cout << "Doing QCD background" << std::endl;
  TH1F* h1 = new TH1F("h1","",1,-10,10);
  tWJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbinSSEnrichW);
  float WsbinSSEnrichW  = h1->Integral()*lumiFact;
  h1->Reset();
  tWJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbinSS);
  float WsbinSS       = h1->Integral()*lumiFact;
  h1->Reset();
  tTTJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbinSSEnrichW);
  float TTsbinSSEnrichW = h1->Integral()*lumiFact;
  h1->Reset();
  tData->Draw("etaL1>>h1",   sbinSSEnrichW);
  float DatasbinSSEnrichW = h1->Integral();
  h1->Reset();
  tData->Draw("etaL1>>h1",   sbinSS);
  float DatasbinSS = h1->Integral();
  h1->Reset();
  
  std::cout << "W events in binSS (MC) = " << WsbinSS << endl; 
  float WsbinSS   = (DatasbinSSEnrichW - TTsbinSSEnrichW)*(WsbinSS/WsbinSSEnrichW);
  std::cout << "W events in sbinSS (corrected) = (" << DatasbinSSEnrichW << " - " << TTsbinSSEnrichW << " )*" << WsbinSS/WsbinSSEnrichW << " = " << WsbinSS << endl;
  float QCDsbinSS = DatasbinSS - WsbinSS;
  std::cout << "SS QCD events = " << QCDsbinSS << std::endl; 
  float QCDsbin   = QCDsbinSS*OStoSSRatio;
  std::cout << "OS QCD events = " << QCDsbin << std::endl;
  /////////////////////////////////////////////////////////////////////

  if(verbose) std::cout << "Doing plots and cuts" << std::endl;
  // Draw with cuts and weights !!!
  if(verbose) std::cout << "Doing OS data plot" << std::endl;
  tData->Draw(  variable+">>hData",   sbin);
  if(verbose) std::cout << "Doing QCD plot" << std::endl;
  tData->Draw(  variable+">>hQCD",    sbinSSRelIso);
  if(verbose) std::cout << "Doing DY plot" << std::endl;
  tDYJets->Draw(variable+">>hDYJets", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  if(verbose) std::cout << "Doing OS W plot" << std::endl;
  tWJets->Draw( variable+">>hWJets",  "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  if(verbose) std::cout << "Doing SS W plot" << std::endl;
  tWJets->Draw( variable+">>hWJetsSS","puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbinSS);
  if(verbose) std::cout << "Doing TTbar plot" << std::endl;
  tTTJets->Draw(variable+">>hTTJets", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  if(verbose) std::cout << "Doing VBF plot" << std::endl;
  tVBFH130->Draw(variable+">>hVBFH130", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  if(verbose) std::cout << "Doing GluGlu plot" << std::endl;
  tGGFH130->Draw(variable+">>hGGFH130", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight*HqTWeight"*sbin);
  
  // Scale histograms
  if(verbose) std::cout << "Scaling histograms" << std::endl;
  hDYJets->Scale( lumiFact );
  hWJets->Scale(     WsbinCorr/hWJets->Integral());
  hWJetsSS->Scale( WsbinSS/hWJetsSS->Integral());
  hTTJets->Scale( lumiFact );
  hQCD->Add(   hWJetsSS, -1);
  hQCD->Scale( QCDsbin/hQCD->Integral());
  hVBFH130->Scale( lumiFact*10 );
  hGGFH130->Scale( lumiFact*10 );
  hVBFH130->Add(hGGFH130,1.0);

  // Add all together
  if(verbose) std::cout << "Stacking histgrams" << std::endl;
  THStack* aStack = new THStack("aStack","");
  aStack->Add(hTTJets);
  aStack->Add(hQCD);
  aStack->Add(hWJets);
  aStack->Add(hDYJets);
  aStack->Add(hVBFH130);

  hData->Sumw2();
  //hData->Draw("P"); // draw signal+bg MC first
  aStack->Draw("HIST");
  hData->Draw("PSAME");

  std::cout << "Data = " << hData->Integral() << endl;
  std::cout << "DYJets = " << hDYJets->Integral() << endl;
  std::cout << "TT = " << hTTJets->Integral() << endl;
  std::cout << "QCD = " << hQCD->Integral() << endl;
  std::cout << "WJets = " << hWJets->Integral() << endl;

  // Legend
  TLegend* leg = new TLegend(0.52,0.50,0.75,0.87,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  leg->SetHeader("#splitline{CMS Preliminary 2011}{#sqrt{s}=7 TeV, L=4.7 fb^{-1}}");
  leg->AddEntry(hData,"Observed","P");
  leg->AddEntry(hDYJets,"Z#rightarrow#tau#tau","F");
  leg->AddEntry(hWJets,"W+jets","F");
  leg->AddEntry(hQCD,"QCD","F");
  leg->AddEntry(hTTJets,"t#bar{t}","F");
  leg->AddEntry(hVBFH130,"10 X H#rightarrow#tau#tau, m_{H}=130 GeV","F");
  leg->Draw();

}
