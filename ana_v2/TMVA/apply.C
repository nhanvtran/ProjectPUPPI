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

void apply(std::string iName="train/OutputTmp.root") { 
  TMVA::Tools::Instance();
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

  float lPt        = 0; reader->AddVariable("pt"                 , &lPt);
  //float lEta       = 0; reader->AddVariable("eta"                , &lEta);
  //float lDR        = 0; reader->AddVariable("dR"                 , &lDR);
  //float lPtc       = 0; reader->AddVariable("ptc"                , &lPtc);
  // float lPtdR      = 0; reader->AddVariable("ptdR"               , &lPtdR);
  //float lPuppi     = 0; reader->AddVariable("puppi"              , &lPuppi);
  float lPtODR     = 0; reader->AddVariable("ptodR"              , &lPtODR);
  //float lPtODRS    = 0; reader->AddVariable("ptodRS"             , &lPtODRS);
  float lPtODRSO   = 0; reader->AddVariable("ptodRSO"            , &lPtODRSO);
  //float lDRLV      = 0; reader->AddVariable("dR_lv"              , &lDRLV);
  //float lPtcLV     = 0; reader->AddVariable("ptc_lv"             , &lPtcLV);
  //float lPtdRLV    = 0; reader->AddVariable("ptdR_lv"            , &lPtdRLV);
  //float lPuppiLV   = 0; reader->AddVariable("puppi_lv"           , &lPuppiLV);
  float lPtODRLV   = 0; reader->AddVariable("ptodR_lv"           , &lPtODRLV);
  //float lPtODRSLV  = 0; reader->AddVariable("ptodRS_lv"          , &lPtODRSLV);
  float lPtODRSOLV = 0; reader->AddVariable("ptodRSO_lv"         , &lPtODRSOLV);
  //float lDRPU      = 0; reader->AddVariable("dR_pu"              , &lDRPU);
  //float lPtcPU     = 0; reader->AddVariable("pt_pu"              , &lPtcPU);
  //float lPtdRPU    = 0; reader->AddVariable("ptdR_pu"            , &lPtdRPU);
  //float lPuppiPU   = 0; reader->AddVariable("puppi_pu"           , &lPuppiPU);
  //float lPtODRPU   = 0; reader->AddVariable("ptodR_pu"           , &lPtODRPU);
  //float lPtODRSPU  = 0; reader->AddVariable("ptodRS_pu"          , &lPtODRSPU);
  //float lPtODRSOPU = 0; reader->AddVariable("ptodRSO_pu"         , &lPtODRSOPU);
  
  std::string lJetName = "BDT";
  reader->BookMVA(lJetName .c_str(),(std::string("weights/TMVAClassificationCategory_PUDisc_v1")+std::string(".weights.xml")).c_str());
  
  TFile *lFile = new TFile(iName.c_str());
  TTree *lTree = (TTree*) lFile->Get("tree");
   lTree->SetBranchAddress("pt"                 , &lPt);
  //lTree->SetBranchAddress("eta"                , &lEta);
   //lTree->SetBranchAddress("dR"                 , &lDR);
   //lTree->SetBranchAddress("ptc"                , &lPtc);
   //lTree->SetBranchAddress("ptdR"               , &lPtdR);
  //lTree->SetBranchAddress("puppi"              , &lPuppi);
  lTree->SetBranchAddress("ptodR"              , &lPtODR);
  //lTree->SetBranchAddress("ptodRS"             , &lPtODRS);
  lTree->SetBranchAddress("ptodRSO"            , &lPtODRSO);
  //lTree->SetBranchAddress("dR_lv"              , &lDRLV);
  // lTree->SetBranchAddress("ptc_lv"             , &lPtcLV);
  //lTree->SetBranchAddress("ptdR_lv"            , &lPtdRLV);
  //lTree->SetBranchAddress("puppi_lv"           , &lPuppiLV);
  lTree->SetBranchAddress("ptodR_lv"           , &lPtODRLV);
  //lTree->SetBranchAddress("ptodRS_lv"          , &lPtODRSLV);
  lTree->SetBranchAddress("ptodRSO_lv"         , &lPtODRSOLV);
  //lTree->SetBranchAddress("dR_pu"              , &lDRPU);
  //lTree->SetBranchAddress("pt_pu"              , &lPtcPU);
  //lTree->SetBranchAddress("ptdR_pu"            , &lPtdRPU);
  //lTree->SetBranchAddress("puppi_pu"           , &lPuppiPU);
  //lTree->SetBranchAddress("ptodR_pu"           , &lPtODRPU);
  //lTree->SetBranchAddress("ptodRS_pu"          , &lPtODRSPU);
  //lTree->SetBranchAddress("ptodRSO_pu"         , &lPtODRSOPU);
    
  int lNEvents = lTree->GetEntries();
  TFile *lOFile = new TFile("Output.root","RECREATE");
  TTree *lOTree = lTree->CloneTree(0);
  float lMVA    = 0; lOTree->Branch("bdt"     ,&lMVA ,"lMVA/F");
  for (Long64_t i0=0; i0<lNEvents;i0++) {
    if (i0 % 10000 == 0) std::cout << "--- ... Processing event: " << double(i0)/double(lNEvents) << std::endl;
    lTree->GetEntry(i0);
    lMVA      = float(reader->EvaluateMVA(lJetName.c_str()));
    lOTree->Fill();
  }
  lOTree->Write();
  lOFile->Close();
  delete reader;
}
