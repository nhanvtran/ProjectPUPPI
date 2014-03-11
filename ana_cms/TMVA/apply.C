#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/Tools.h"
#endif

void apply(std::string iName="../RecoOutput_vLoose.root") { //train/OutputTmp.root") { 
  TMVA::Tools::Instance();
  TMVA::Reader *reader  = new TMVA::Reader( "!Color:!Silent" );    
  TMVA::Reader *reader2 = new TMVA::Reader( "!Color:!Silent" );    
 
  float lTK        = 0; reader->AddVariable("trk"                , &lTK);
  float lVtx       = 0; reader->AddVariable("vtx"                , &lVtx);
  float lVId       = 0; reader->AddVariable("vtxid"              , &lVId);
  float lD0        = 0; reader->AddVariable("d0"                 , &lD0);
  float lDZ        = 0; reader->AddVariable("dZ"                 , &lDZ);
  float lDepth     = 0; reader->AddVariable("depth"              , &lDepth);
  float lTime      = 0; reader->AddVariable("time"               , &lTime);
  reader2->AddVariable("depth"              , &lDepth);
  reader2->AddVariable("time"               , &lTime);
  
  std::string lJetName = "BDT";
  reader ->BookMVA("type1",(std::string("weights/TMVAClassificationCategory_PUDisc_exp_pftype1")+std::string(".weights.xml")).c_str());
  reader2->BookMVA("type4",(std::string("weights/TMVAClassificationCategory_PUDisc_exp_pftype4")+std::string(".weights.xml")).c_str());
  reader2->BookMVA("type5",(std::string("weights/TMVAClassificationCategory_PUDisc_exp_pftype5")+std::string(".weights.xml")).c_str());
  
  TFile *lFile = new TFile(iName.c_str());
  TTree *lTree = (TTree*) lFile->Get("Result");
  lTree->SetBranchAddress("trk"                , &lTK);
  lTree->SetBranchAddress("vtx"                , &lVtx);
  int lIVId = 0; lTree->SetBranchAddress("vtxid"               , &lIVId);
  lTree->SetBranchAddress("d0"                , &lD0);
  lTree->SetBranchAddress("dZ"                , &lDZ);
  lTree->SetBranchAddress("depth"             , &lDepth);
  lTree->SetBranchAddress("time"              , &lTime);
  int lPFType = 0; lTree->SetBranchAddress("pftype"            , &lPFType);
    
  int lNEvents = lTree->GetEntries();
  TFile *lOFile = new TFile("Output.root","RECREATE");
  TTree *lOTree = lTree->CloneTree(0);
  float lMVA    = 0; lOTree->Branch("bdt"     ,&lMVA ,"lMVA/F");
  for (Long64_t i0=0; i0<lNEvents;i0++) {
    if (i0 % 10000 == 0) std::cout << "--- ... Processing event: " << double(i0)/double(lNEvents) << std::endl;
    lTree->GetEntry(i0);
    lVId = float(lIVId);
    lMVA = 1;
    if(lPFType == 1 || lPFType == 2 || lPFType == 3) lMVA      = float(reader->EvaluateMVA("type1"));
    if(lPFType == 4) lMVA      = float(reader2->EvaluateMVA("type4"));
    if(lPFType == 5) lMVA      = float(reader2->EvaluateMVA("type5"));
    lOTree->Fill();
  }
  lOTree->Write();
  lOFile->Close();
  delete reader;
}
