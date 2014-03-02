#include "GenMatcher.hh"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <iostream>

GenMatcher *fMatch = 0; 

TTree* load(std::string iName) { 
  TFile *lFile = new TFile(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  return lTree;
}
int main( int argc, char **argv ) {
  gROOT->ProcessLine("#include <vector>");          
  bool         lGen = atoi(argv[3]);
  std::cout << "===> Gen " << lGen << std::endl;
  std::string lName = argv[2];
  int maxEvents     = atoi(argv[1]);
  TTree *lTree = load(lName); 
  if(lTree->GetEntries() < maxEvents || maxEvents == -1) maxEvents = lTree->GetEntries(); 
  fMatch    = new GenMatcher(lTree);
  TFile *lFile = new TFile("Output.root","RECREATE");
  TTree *lOut  = new TTree("Result","Result");
  fMatch->setup(lOut);
  for(int i0 = 0; i0 < maxEvents; i0++) { 
    if(i0 % 50 == 0) std::cout << "===> Processed " << i0 << " - Done : " << (float(i0)/float(maxEvents)) << std::endl;
    fMatch   ->load(i0);
    if(!lGen) fMatch->match();
    if(lGen ) fMatch->fillGen();
  }
  lFile->cd();
  lOut->Write();
  lFile->Close();
}
