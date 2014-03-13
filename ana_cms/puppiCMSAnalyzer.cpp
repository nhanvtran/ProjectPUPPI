#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include "TROOT.h"
#include "TFile.h"
#include "TList.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TEllipse.h"
#include "Math/ProbFunc.h"

#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

#include "puppiContainer.hh"
#include "NoTrees.hh"
#include "RecoObj.hh"

using namespace std;
using namespace fastjet;

/////////////////////////////////////////////////////////////////////
// Tree variables
/////////////////////////////////////////////////////////////////////
int   fCount      = 0; 
float fPartPt     = 0; 
float fPartEta    = 0; 
float fPartPhi    = 0; 
float fPartM      = 0;
float fPartTk     = 0; 
float fPartVtx    = 0; 
int   fPartPFType = 0; 
float fPartDepth  = 0; 
float fPartTime   = 0; 
float fGenPt      = 0; 
int   fVtxId      = 0; 

float fPt     = 0; 
float fEta    = 0; 
float fPhi    = 0; 
float fM      = 0;
float fDPhi   = 0;
float fDEta   = 0; 
float fDR     = 0; 
float fCh     = 0;
float fPU     = 0;
float fVtx    = 0; 
float fTk     = 0;
int   fPFType = 0; 
float fDepth  = 0; 
float fTime   = 0; 
int   fVId    = 0;  
float fW0     = 0; 
float fW1     = 0; 
float fW2     = 0; 
float fW3     = 0; 
float fW4     = 0; 
float fW5     = 0; 
float fW6     = 0; 

float fW0L     = 0; 
float fW1L     = 0; 
float fW2L     = 0; 
float fW3L     = 0; 
float fW4L     = 0; 
float fW5L     = 0; 
float fW6L     = 0; 

float fW0P  = 0; 
float fW1P  = 0; 
float fW2P  = 0; 
float fW3P  = 0; 
float fW4P  = 0; 
float fW5P  = 0; 
float fW6P  = 0; 

TTree *fTree;

/////////////////////////////////////////////////////////////////////
// Helper functions
/////////////////////////////////////////////////////////////////////
void setupCMSReadOut(TTree *iTree );
void readCMSEvent( TTree *iTree, std::vector< RecoObj > &allParticles);
void readEvent( std::vector< fastjet::PseudoJet > &allParticles, std::vector<int> &v_isPU, std::vector<int> &v_isCh );
void getJets(std::vector < fastjet::PseudoJet > constits,std::vector < fastjet::PseudoJet > &jets);
void addBranches( TTree &iTree){
  iTree.Branch("pt"         ,&fPt ,"fPt/F");
  iTree.Branch("eta"        ,&fEta,"fEta/F");
  iTree.Branch("phi"        ,&fPhi,"fPhi/F");
  iTree.Branch("m"          ,&fM  ,"fM/F");
  iTree.Branch("dEta"       ,&fDEta   ,"fDEta/F");
  iTree.Branch("dPhi"       ,&fDPhi   ,"fDPhi/F");
  iTree.Branch("dr"         ,&fDR     ,"fDR/F");
  iTree.Branch("charge"     ,&fCh     ,"fCh/F");
  iTree.Branch("pu"         ,&fPU     ,"fPU/F");
  iTree.Branch("time"       ,&fTime   ,"fTime/F");
  iTree.Branch("pfType"     ,&fPFType ,"fPFType/I");
  iTree.Branch("depth"      ,&fDepth  ,"fDepth/F");
  iTree.Branch("tk"         ,&fTk ,"fTk/F");
  iTree.Branch("vtx"        ,&fVtx ,"fVtx/F");
  iTree.Branch("vid"        ,&fVId ,"fVId/I");
  iTree.Branch("dR"         ,&fW0 ,"fW0/F");								    
  iTree.Branch("ptc"        ,&fW1 ,"fW1/F");								    
  iTree.Branch("ptdR"       ,&fW2 ,"fW2/F");								    
  iTree.Branch("puppi"      ,&fW3 ,"fW3/F");								    
  iTree.Branch("ptodR"      ,&fW4 ,"fW4/F");								    
  iTree.Branch("ptodRS"     ,&fW5 ,"fW5/F");								    
  iTree.Branch("ptodRSO"    ,&fW6 ,"fW6/F");								    

  iTree.Branch("dR_lv"         ,&fW0L ,"fW0L/F");								    
  iTree.Branch("ptc_lv"        ,&fW1L ,"fW1L/F");								    
  iTree.Branch("ptdR_lv"       ,&fW2L ,"fW2L/F");								    
  iTree.Branch("puppi_lv"      ,&fW3L ,"fW3L/F");								    
  iTree.Branch("ptodR_lv"      ,&fW4L ,"fW4L/F");								    
  iTree.Branch("ptodRS_lv"     ,&fW5L ,"fW5L/F");								    
  iTree.Branch("ptodRSO_lv"    ,&fW6L ,"fW6L/F");								    

  iTree.Branch("dR_pu"         ,&fW0P ,"fW0P/F");								    
  iTree.Branch("ptc_pu"        ,&fW1P ,"fW1P/F");								    
  iTree.Branch("ptdR_pu"       ,&fW2P ,"fW2P/F");								    
  iTree.Branch("puppi_pu"      ,&fW3P ,"fW3P/F");								    
  iTree.Branch("ptodR_pu"      ,&fW4P ,"fW4P/F");								    
  iTree.Branch("ptodRS_pu"     ,&fW5P ,"fW5P/F");								    
  iTree.Branch("ptodRSO_pu"    ,&fW6P ,"fW6P/F");								    
}
int main( int argc, char **argv ) {
  bool iCorr = false;
  std::string lName = argv[2];
  int maxEvents     = atoi(argv[1]);
  gROOT->ProcessLine("#include <vector>");                
  int nEvts = 0;
  //int maxEvents = 1000;
  
  TFile* lReadout = new TFile("RecoOutput_vMedium.root");
  TTree* fTree    = (TTree*) lReadout->FindObjectAny("Result");
  setupCMSReadOut(fTree);
  TFile* lFile = new TFile("output/outputTmp.root", "RECREATE");
  TTree* lTree      = new TTree("tree", "tree");
  addBranches(*lTree);
  std::vector < RecoObj > allParticles;
  lReadout->cd();
  readCMSEvent(fTree, allParticles);        
  while(true){
      nEvts++;
      if (nEvts % 10 == 0) std::cout << "event no. = " << nEvts << std::endl;
      if (nEvts == maxEvents){ break; }
      lReadout->cd();
      readCMSEvent(fTree, allParticles );        
      puppiContainer curEvent(allParticles);
      std::vector<fastjet::PseudoJet> puParticles = curEvent.puParticles();
      std::vector<fastjet::PseudoJet> pvParticles = curEvent.pvParticles();
      std::vector<fastjet::PseudoJet> pfParticles = curEvent.pfParticles();
      for(unsigned int i0 = 0; i0 < allParticles.size(); i0++ ) { 
	fPt      = allParticles[i0].pt; 
	fEta    = allParticles[i0].eta; 
	fPhi    = allParticles[i0].phi; 
	//fDR   = sqrt(fDPhi*fDPhi + fDEta*fDEta);
	fPU     = allParticles[i0].id > 3;
	fCh     = (allParticles[i0].id > 1 && allParticles[i0].id < 4);
	fTk     = allParticles[i0].trkChi2;
	fVtx    = allParticles[i0].vtxChi2;
	fVId    = allParticles[i0].vtxId;
	fPFType = allParticles[i0].pfType;
	fDepth  = allParticles[i0].depth;
	fTime   = allParticles[i0].time;
	//Should make the bottom a vector, but I'm lazy
	fW0   = curEvent.goodVar(i0,pfParticles,0);//compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,0),lAvgPU[0] ,lAvgPU[1]);
	fW1   = curEvent.goodVar(i0,pfParticles,1);//compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,1),lAvgPU[2] ,lAvgPU[3]);
	fW2   = curEvent.goodVar(i0,pfParticles,2);//compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,2),lAvgPU[4] ,lAvgPU[5]);
	fW3   = curEvent.goodVar(i0,pfParticles,7);//compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,3),lAvgPU[6] ,lAvgPU[7]);
	fW4   = curEvent.goodVar(i0,pfParticles,4);// compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,4),lAvgPU[8] ,lAvgPU[9]);
	fW5   = curEvent.goodVar(i0,pfParticles,8);// compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,5),lAvgPU[10],lAvgPU[11]);
	fW6   = curEvent.goodVar(i0,pfParticles,6);// compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,6),lAvgPU[12],lAvgPU[13]);

	fW0L  = curEvent.goodVar(i0,pvParticles,0);//compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,0),lAvgPU[14],lAvgPU[15]);
	fW1L  = curEvent.goodVar(i0,pvParticles,1);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,1),lAvgPU[16],lAvgPU[17]);
	fW2L  = curEvent.goodVar(i0,pvParticles,2);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,2),lAvgPU[18],lAvgPU[19]);
	fW3L  = curEvent.goodVar(i0,pvParticles,7);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,3),lAvgPU[20],lAvgPU[21]);
	fW4L  = curEvent.goodVar(i0,pvParticles,4);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,4),lAvgPU[22],lAvgPU[23]);
	fW5L  = curEvent.goodVar(i0,pvParticles,8);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,5),lAvgPU[24],lAvgPU[25]);
	fW6L  = curEvent.goodVar(i0,pvParticles,6);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,6),lAvgPU[26],lAvgPU[27]);

	fW0P  = curEvent.goodVar(i0,puParticles,0);//compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,0),lAvgPU[14],lAvgPU[15]);
	fW1P  = curEvent.goodVar(i0,puParticles,1);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,1),lAvgPU[16],lAvgPU[17]);
	fW2P  = curEvent.goodVar(i0,puParticles,2);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,2),lAvgPU[18],lAvgPU[19]);
	fW3P  = curEvent.goodVar(i0,puParticles,7);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,3),lAvgPU[20],lAvgPU[21]);
	fW4P  = curEvent.goodVar(i0,puParticles,4);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,4),lAvgPU[22],lAvgPU[23]);
	fW5P  = curEvent.goodVar(i0,puParticles,8);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,5),lAvgPU[24],lAvgPU[25]);
	fW6P  = curEvent.goodVar(i0,puParticles,6);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,6),lAvgPU[26],lAvgPU[27]);
	lTree->Fill();
      }
      allParticles.clear();
  }
  lFile->cd();
  lTree->Write();
}
void setupCMSReadOut(TTree *iTree ) { 
  fCount = 0;
  iTree->SetBranchAddress("pt"    ,&fPartPt);
  iTree->SetBranchAddress("eta"   ,&fPartEta);
  iTree->SetBranchAddress("phi"   ,&fPartPhi);
  iTree->SetBranchAddress("m"     ,&fPartM);
  iTree->SetBranchAddress("genpt" ,&fGenPt);
  iTree->SetBranchAddress("vtxid" ,&fVtxId);
  iTree->SetBranchAddress("vtx"   ,&fPartVtx);
  iTree->SetBranchAddress("trk"   ,&fPartTk);
  iTree->SetBranchAddress("pftype",&fPartPFType);  
  iTree->SetBranchAddress("depth" ,&fPartDepth);
  iTree->SetBranchAddress("time"  ,&fPartTime);
}
void readCMSEvent(TTree *iTree, std::vector< RecoObj > &allParticles) { 
  float px, py, pz, e, pdgid, isCh, isPU,isPV = 0;
  while(true){
    iTree->GetEntry(fCount);
    fCount++;
    if(fPartPt == -1) break;
    isCh    = (fPartPFType == 1 || fPartPFType == 2 || fPartPFType == 3);
    isPV   = (fVtxId == 0);
    isPU   = (fGenPt/fPartPt < 0.);
    //if(fGenPt/fPartPt > 0. && fGenPt/fPartPt < 0.8) continue;
    TLorentzVector pVec; pVec.SetPtEtaPhiM(fPartPt,fPartEta,fPartPhi,fPartM);
    int lId = 0; if(!isPV) lId++; if(isCh) lId+=2;
    if(isPU) lId += 4;
    RecoObj curPseudoJet;
    curPseudoJet.pt  = fPartPt;
    curPseudoJet.eta = fPartEta;
    curPseudoJet.phi = fPartPhi;
    curPseudoJet.m   = fPartM;
    curPseudoJet.id  = lId;
    curPseudoJet.vtxId = fVtxId;
    curPseudoJet.trkChi2 = fPartTk;
    curPseudoJet.vtxChi2 = fPartVtx;
    curPseudoJet.pfType  = fPartPFType;
    curPseudoJet.depth   = fPartDepth;
    curPseudoJet.time    = fPartTime;
    allParticles.push_back( curPseudoJet );
  }
}
void getJets(std::vector < fastjet::PseudoJet > constits,std::vector < fastjet::PseudoJet > &jets) { 
  double rParam = 0.7;
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);    
  int activeAreaRepeats = 1;
  double ghostArea = 0.01;
  double ghostEtaMax = 7.0;
  fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
  fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
  fastjet::ClusterSequenceArea* thisClustering_ = new fastjet::ClusterSequenceArea(constits, jetDef, fjAreaDefinition);
  std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering_->inclusive_jets(5.0));
  for(unsigned int i0 = 0; i0 < out_jets.size(); i0++) jets.push_back(out_jets[i0]);
}
