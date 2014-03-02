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
#include "TMVA/Factory.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/Tools.h"
#endif

void classify() {
  TMVA::Tools::Instance();
  std::cout << "==> Start TMVAClassification" << std::endl;

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassificationCategory", outputFile,
					       "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   // A very simple MVA (feel free to uncomment and comment what you like) => as a rule of thumb 10-20 variables is where people start to get worried about total number

   factory->AddVariable("pt"         ,'F');
   factory->AddVariable("eta"        ,'F');
   factory->AddVariable("tk"         ,'F');
   factory->AddVariable("vtx"        ,'F');
   factory->AddVariable("vid"        ,'I');
   //factory->AddVariable("dR"         ,'F');
   //factory->AddVariable("ptc"        ,'F');
   //factory->AddVariable("ptdR"       ,'F');
   //factory->AddVariable("puppi"      ,'F');
   factory->AddVariable("ptodR"      ,'F');
   //factory->AddVariable("ptodRS"     ,'F');
   factory->AddVariable("ptodRSO"    ,'F');

   //factory->AddVariable("dR_lv"         ,'F');
   //factory->AddVariable("ptc_lv"        ,'F');
   //factory->AddVariable("ptdR_lv"       ,'F');
   //factory->AddVariable("puppi_lv"      ,'F');
   factory->AddVariable("ptodR_lv"      ,'F');
   //factory->AddVariable("ptodRS_lv"     ,'F');
   factory->AddVariable("ptodRSO_lv"    ,'F');

   //factory->AddVariable("dR_pu"         ,'F');
   //factory->AddVariable("pt_pu"         ,'F');
   //factory->AddVariable("ptdR_pu"       ,'F');
   //factory->AddVariable("puppi_pu"      ,'F');
   //factory->AddVariable("ptodR_pu"      ,'F');
   //factory->AddVariable("ptodRS_pu"     ,'F');
   //factory->AddVariable("ptodRSO_pu"    ,'F');
   
   TCut lCut   = ""; //No Cut needed

   TString lBase         = "train/";
   TString lSAName       = lBase; lSAName+="signal.root";  TFile *lSAInput = TFile::Open( lSAName );
   TTree   *lSASignal    = (TTree*)lSAInput    ->Get("tree"); 
   TString lSBName       = lBase; lSBName+="pu.root";    TFile *lSBInput = TFile::Open( lSBName );
   TTree   *lSBSignal    = (TTree*)lSBInput    ->Get("tree"); 
   
   Double_t lSWeight = 1.0;
   Double_t lBWeight = 1.0;
   gROOT->cd( outfileName+TString(":/") );   
   factory->AddSignalTree    ( lSASignal, lSWeight );
   
   gROOT->cd( outfileName+TString(":/") );   
   factory->AddBackgroundTree( lSBSignal, lBWeight );
   
   factory->PrepareTrainingAndTestTree( lCut, lCut,"nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Block:NormMode=NumEvents:!V" );
   TString lBDTDef   = "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=F:nCuts=2000:NNodesMax=10000:MaxDepth=20:SeparationType=GiniIndex";
   //TString lBDTDef   = "!H:!V:NTrees=100:BoostType=Grad:Shrinkage=0.05:UseBaggedGrad=F:nCuts=2000:NNodesMax=10000:MaxDepth=5:UseYesNoLeaf=F:nEventsMin=200";
   factory->BookMethod(TMVA::Types::kBDT,"PUDisc_v1",lBDTDef);   
   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;
   delete factory;
   //if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
