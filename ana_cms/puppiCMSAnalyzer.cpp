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

float fW0P     = 0; 
float fW1P     = 0;
float fW2P     = 0; 
float fW3P     = 0; 
float fW4P     = 0; 
float fW5P     = 0; 
float fW6P     = 0; 

TTree *fTree;

/////////////////////////////////////////////////////////////////////
// Helper functions
/////////////////////////////////////////////////////////////////////
void setupCMSReadOut(TTree *iTree );
void readCMSEvent   ( TTree *iTree, std::vector< fastjet::PseudoJet > &allParticles, 
		      std::vector<int> &v_isPU, std::vector<int> &v_isCh,std::vector<float> &v_pTk,std::vector<float> &v_pVtx,std::vector<int> &v_pVId );
void readEvent( std::vector< fastjet::PseudoJet > &allParticles, std::vector<int> &v_isPU, std::vector<int> &v_isCh );
void getJets  (std::vector < fastjet::PseudoJet > constits,std::vector < fastjet::PseudoJet > &jets);
void addBranches( TTree &iTree){
  iTree.Branch("pt"         ,&fPt ,"fPt/F");
  iTree.Branch("eta"        ,&fEta,"fEta/F");
  iTree.Branch("phi"        ,&fPhi,"fPhi/F");
  iTree.Branch("m"          ,&fM  ,"fM/F");
  iTree.Branch("dEta"       ,&fDEta ,"fDEta/F");
  iTree.Branch("dPhi"       ,&fDPhi ,"fDPhi/F");
  iTree.Branch("dr"         ,&fDR ,"fDR/F");
  iTree.Branch("charge"     ,&fCh ,"fCh/F");
  iTree.Branch("pu"         ,&fPU ,"fPU/F");
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
  iTree.Branch("pt_pu"         ,&fW1P ,"fW1P/F");								    
  iTree.Branch("ptdR_pu"       ,&fW2P ,"fW2P/F");								    
  iTree.Branch("puppi_pu"      ,&fW3P ,"fW3P/F");								    
  iTree.Branch("ptodR_pu"      ,&fW4P ,"fW4P/F");								    
  iTree.Branch("ptodRS_pu"     ,&fW5P ,"fW5P/F");								    
  iTree.Branch("ptodRSO_pu"    ,&fW6P ,"fW6P/F");								    
}
double* getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,std::vector<int> &iIsPU,puppiContainer &iEvent) { 
  std::vector<double> lValsPV;
  std::vector<double> lValsPU;
  for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) { 
    if( iConstits[i0].pt() < 0.5) continue;
    double pVal = iEvent.goodVar(iConstits[i0],iParticles,iOpt);
    if(iConstits[i0].user_index() % 2 == 1) lValsPU.push_back(pVal);
    if(iConstits[i0].user_index() % 2 == 0) lValsPV.push_back(pVal);
    //if( iIsPU[iConstits[i0].user_index()]) lAvgPU += pVal;
  }
  std::sort (lValsPV.begin(),lValsPV.end());   
  std::sort (lValsPU.begin(),lValsPU.end());   
  double lMedPV = 0; if(lValsPV.size() > 0) lMedPV = lValsPV[int(lValsPV.size()*0.5+0.5)];
  double lMedPU = 0; if(lValsPU.size() > 0) lMedPU = lValsPU[int(lValsPU.size()*0.5+0.5)];
  double lRMSPU = 0; 
 //if(iOpt == 6) cout << "===> Med : " << lMedPU << " -- " << lAvgPU/lValsPU.size() << endl;
  //lMedPU = lAvgPU/lValsPU.size();
  for(unsigned int i0 = 0 ;i0 < lValsPU.size(); i0++) {lRMSPU += (lValsPU[i0]-lMedPU)*(lValsPU[i0]-lMedPU);}
  double lRMSPV = 0; 
  for(unsigned int i0 = 0 ;i0 < lValsPV.size(); i0++) {lRMSPV += (lValsPV[i0]-lMedPV)*(lValsPV[i0]-lMedPV);}
  if(lValsPV.size() > 0)  lRMSPV/=lValsPV.size(); 
  if(lValsPU.size() > 0)  lRMSPU/=lValsPU.size(); 
   double *lVals = new double[4];
  lVals[0] = lMedPV;  
  lVals[1] = sqrt(lRMSPV);
  lVals[2] = lMedPU;  
  lVals[3] = sqrt(lRMSPU);
  return lVals;
}
void getAvgVector(std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,std::vector<int> &iIsPU,puppiContainer &iEvent,
		  std::vector<double> &iAvgPV,std::vector<double> &iAvgPU ) { 
  //Compute the average and RMS store thems sequentially
  for(int i0 = 0; i0 < 7; i0++) { 
    double *pVals = getRMSAvg(i0,iConstits,iParticles,iIsPU,iEvent); 
    for(int i1 = 0; i1 < 2; i1++) {
      iAvgPV.push_back(pVals[i1]); 
      iAvgPU.push_back(pVals[i1+2]);
    }
  }
}
double compute(bool iCorr,double iVal,double iMed,double iRMS) { 
  if(!iCorr) return iVal;
  double lVal = (iVal-iMed)/iRMS;
  return lVal;
}
int main( int argc, char **argv ) {
  bool iCorr = false;
  std::string lName = argv[2];
  int maxEvents     = atoi(argv[1]);
  gROOT->ProcessLine("#include <vector>");                
  int nEvts = 0;
  
  TFile* lReadout = new TFile("Output.root");
  TTree* fTree    = (TTree*) lReadout->FindObjectAny("Result");
  setupCMSReadOut(fTree);
  TFile* lFile = new TFile("output/outputTmp.root", "RECREATE");
  TTree* lTree      = new TTree("tree", "tree");
  addBranches(*lTree);
  std::vector < fastjet::PseudoJet > allParticles;
  std::vector < int   > v_isPU;
  std::vector < int   > v_isCh;
  std::vector < float > v_pTk;
  std::vector < float > v_pVtx;
  std::vector < int   > v_pVId;
  lReadout->cd();
  readCMSEvent(fTree, allParticles, v_isPU, v_isCh,v_pTk,v_pVtx,v_pVId );        
  while(true){
      nEvts++;
      if (nEvts % 10 == 0) std::cout << "event no. = " << nEvts << std::endl;
      if (nEvts == maxEvents){ break; }
      lReadout->cd();
      readCMSEvent(fTree, allParticles, v_isPU, v_isCh,v_pTk,v_pVtx,v_pVId );        
      puppiContainer curEvent(allParticles, v_isPU, v_isCh,0);
      std::vector<fastjet::PseudoJet> pvParticles;
      for(int i0 = 0; i0 < allParticles.size(); i0++) {
	if(v_pVId[i0] == 0) pvParticles.push_back(allParticles[i0]);
      }
      // getJets(curEvent.pfParticles(),jets);
      //std::vector<PseudoJet> lConstits = curEvent.pfParticles();
      //std::vector<double>  lAvgPV ; //store the Average and RMS for PV
      //std::vector<double>  lAvgPU ; //Ditto for the PU 
      //if(iCorr) { 
	//fastjet::Selector lSel = fastjet::SelectorCircle(1.0);
	//lSel.set_reference(jets[0]);
	//std::vector<PseudoJet> lNearParticles = lSel(curEvent.pfParticles());
	//getAvgVector(lNearParticles  ,allParticles,v_isPU,curEvent,lAvgPV,lAvgPU);
	//getAvgVector(allParticles  ,pvParticles ,v_isPU,curEvent,lAvgPV,lAvgPU);
	//getAvgVector(lNearParticles,puParticles ,v_isPU,curEvent,lAvgPV,lAvgPU);
	// }
      //lFile->cd();
      for(unsigned int i0 = 0; i0 < allParticles.size(); i0++ ) { 
	if( allParticles[i0].pt() < 0.5) continue;
	fPt   = allParticles[i0].pt(); 
	fEta  = allParticles[i0].eta(); 
	//fDPhi = (jets[0].phi()-allParticles[i0].phi());
	//if(fabs(fDPhi) > 2.*TMath::Pi()-fabs(fDPhi) && fDPhi < 0.) fDPhi =   2.*TMath::Pi() + fDPhi;
	//if(fabs(fDPhi) > 2.*TMath::Pi()-fabs(fDPhi) && fDPhi > 0.) fDPhi =  -2.*TMath::Pi() + fDPhi;
	//fDEta = (jets[0].eta()-lConstits[i0].eta());
	fDR   = sqrt(fDPhi*fDPhi + fDEta*fDEta);
	fPU   = allParticles[i0].user_index() % 2;
	fCh   = (allParticles[i0].user_index() > 1 && allParticles[i0].user_index() < 4);
	fTk   = v_pTk [i0];
	fVtx  = v_pVtx[i0];
	fVId  = v_pVId[i0];

	//Should make the bottom a vector, but I'm lazy
	fW0   = curEvent.goodVar(allParticles[i0],allParticles,0);//compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,0),lAvgPU[0] ,lAvgPU[1]);
	fW1   = curEvent.goodVar(allParticles[i0],allParticles,1);//compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,1),lAvgPU[2] ,lAvgPU[3]);
	fW2   = curEvent.goodVar(allParticles[i0],allParticles,2);//compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,2),lAvgPU[4] ,lAvgPU[5]);
	fW3   = curEvent.goodVar(allParticles[i0],allParticles,3);//compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,3),lAvgPU[6] ,lAvgPU[7]);
	fW4   = curEvent.goodVar(allParticles[i0],allParticles,4);// compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,4),lAvgPU[8] ,lAvgPU[9]);
	fW5   = curEvent.goodVar(allParticles[i0],allParticles,5);// compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,5),lAvgPU[10],lAvgPU[11]);
	fW6   = curEvent.goodVar(allParticles[i0],allParticles,6);// compute(iCorr,curEvent.goodVar(lConstits[i0],allParticles,6),lAvgPU[12],lAvgPU[13]);

	fW0L  = curEvent.goodVar(allParticles[i0],pvParticles,0);//compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,0),lAvgPU[14],lAvgPU[15]);
	fW1L  = curEvent.goodVar(allParticles[i0],pvParticles,1);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,1),lAvgPU[16],lAvgPU[17]);
	fW2L  = curEvent.goodVar(allParticles[i0],pvParticles,2);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,2),lAvgPU[18],lAvgPU[19]);
	fW3L  = curEvent.goodVar(allParticles[i0],pvParticles,3);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,3),lAvgPU[20],lAvgPU[21]);
	fW4L  = curEvent.goodVar(allParticles[i0],pvParticles,4);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,4),lAvgPU[22],lAvgPU[23]);
	fW5L  = curEvent.goodVar(allParticles[i0],pvParticles,5);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,5),lAvgPU[24],lAvgPU[25]);
	fW6L  = curEvent.goodVar(allParticles[i0],pvParticles,6);// compute(iCorr,curEvent.goodVar(lConstits[i0],pvParticles,6),lAvgPU[26],lAvgPU[27]);
	lTree->Fill();
	continue;
	/*
	fW0P  = compute(iCorr,curEvent.goodVar(lConstits[i0],puParticles,0),lAvgPU[28],lAvgPU[29]);
	fW1P  = compute(iCorr,curEvent.goodVar(lConstits[i0],puParticles,1),lAvgPU[30],lAvgPU[31]);
	fW2P  = compute(iCorr,curEvent.goodVar(lConstits[i0],puParticles,2),lAvgPU[32],lAvgPU[33]);
	fW3P  = compute(iCorr,curEvent.goodVar(lConstits[i0],puParticles,3),lAvgPU[34],lAvgPU[35]);
	fW4P  = compute(iCorr,curEvent.goodVar(lConstits[i0],puParticles,4),lAvgPU[36],lAvgPU[37]);
	fW5P  = compute(iCorr,curEvent.goodVar(lConstits[i0],puParticles,5),lAvgPU[38],lAvgPU[39]);
	fW6P  = compute(iCorr,curEvent.goodVar(lConstits[i0],puParticles,6),lAvgPU[40],lAvgPU[41]);
	*/
      }
      allParticles.clear();
      v_isPU.clear();
      v_isCh.clear();
      v_pTk .clear();
      v_pVtx.clear();
      v_pVId.clear();
  }
  lFile->cd();
  lTree->Write();
}
void setupCMSReadOut(TTree *iTree ) { 
  fCount = 0;
  iTree->SetBranchAddress("pt"   ,&fPartPt);
  iTree->SetBranchAddress("eta"  ,&fPartEta);
  iTree->SetBranchAddress("phi"  ,&fPartPhi);
  iTree->SetBranchAddress("m"    ,&fPartM);
  iTree->SetBranchAddress("genpt",&fGenPt);
  iTree->SetBranchAddress("vtxid",&fVtxId);
  iTree->SetBranchAddress("vtx"  ,&fPartVtx);
  iTree->SetBranchAddress("trk"  ,&fPartTk);
}
void readCMSEvent(TTree *iTree, std::vector< fastjet::PseudoJet > &allParticles, std::vector<int> &v_isPU, std::vector<int> &v_isCh,
		  std::vector<float> &v_pTk,std::vector<float> &v_pVtx,std::vector<int> &v_pVId ){
  float px, py, pz, e, pdgid, isCh, isPU = 0;
  while(true){
    iTree->GetEntry(fCount);
    fCount++;
    if(fPartPt == -1) break;
    isCh = (fVtxId != -1);
    isPU = (fGenPt/fPartPt < 0.8);  //Matching Definition
    TLorentzVector pVec; pVec.SetPtEtaPhiM(fPartPt,fPartEta,fPartPhi,fPartM);
    px = pVec.Px();
    py = pVec.Py();
    pz = pVec.Pz();
    e  = pVec.E();
    if (px == 0 && py == 0 && pz == 0 && e == 0) return;
    // fill vector of pseudojets
    fastjet::PseudoJet curPseudoJet( px, py, pz, e );
    if (fabs(curPseudoJet.eta()) < 5){
      int lId = 0; if(isPU) lId++; if(isCh) lId+=2;
      curPseudoJet.set_user_index(lId);
      allParticles.push_back( curPseudoJet );
      v_isPU.push_back(isPU);
      v_isCh.push_back(isCh);         
      v_pTk .push_back(fPartTk);         
      v_pVtx.push_back(fPartVtx);  
      v_pVId.push_back(fVtxId);
    }
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
