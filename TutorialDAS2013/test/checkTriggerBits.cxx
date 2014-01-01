#include "TTree.h"
#include "TFile.h"
#include <iostream>

void checkTriggerBits(char * inFileName) {

  TFile file(inFileName);
  TTree *tr = (TTree*)file.Get("muTauStreamAnalyzerRaw/tree");

  int nEvents = tr->GetEntries();
  nEvents = 10;
  
  std::vector< int >* tauXTriggers  = new std::vector< int >();
  tr->SetBranchAddress("tauXTriggers",    &tauXTriggers);
  
  for( int iEvent = 0; iEvent < nEvents; iEvent++) {
    tr->GetEntry(iEvent);
    std::cout << "Event " << iEvent << std::endl;
    for( unsigned int iXT = 0; iXT < tauXTriggers->size(); iXT++) {
      std::cout << "[" << iXT << "] " << (*tauXTriggers)[iXT] << std::endl;
    }
  }
}
