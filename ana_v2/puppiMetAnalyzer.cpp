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

#include "fastjet/PseudoJet.hh"

#include "puppiCleanContainer.hh"
#include "puppiTMVAContainer.hh"
#include "NoTrees.hh"


using namespace std;
using namespace fastjet;

ifstream fin;
/////////////////////////////////////////////////////////////////////
// Tree variables
/////////////////////////////////////////////////////////////////////
float fMet      = 0;
float fMetPhi   = 0;
float fSumet    = 0;
float fGMet      = 0;
float fGMetPhi   = 0;
float fGSumet    = 0;
float fTMet      = 0;
float fTMetPhi   = 0;
float fTSumet    = 0;
float fPMet      = 0;
float fPMetPhi   = 0;
float fPSumet    = 0;
/////////////////////////////////////////////////////////////////////
// Helper functions
/////////////////////////////////////////////////////////////////////
void readEvent( std::vector< fastjet::PseudoJet > &allParticles, std::vector<int> &v_isPU, std::vector<int> &v_isCh );
void analyzeEvent(std::vector < fastjet::PseudoJet > &constits, std::vector < fastjet::PseudoJet > &trimconstits, std::vector < fastjet::PseudoJet > &genconstits,std::vector < fastjet::PseudoJet > &puppiconstits);
float *met( std::vector< fastjet::PseudoJet > &iParticles);
void addBranches( TTree &iTree){
  iTree.Branch("genmet"   ,&fGMet    ,"fGMet/F");
  iTree.Branch("genmetphi",&fGMetPhi ,"fGMetPhi/F");
  iTree.Branch("gensumet" ,&fGSumet  ,"fGSumet/F");
  iTree.Branch("met"      ,&fMet     ,"fMet/F");
  iTree.Branch("metphi"   ,&fMetPhi  ,"fMetPhi/F");
  iTree.Branch("sumet"    ,&fSumet   ,"fSumet/F");
  iTree.Branch("trimmet"  ,&fTMet    ,"fTMet/F");
  iTree.Branch("trimmetphi" ,&fTMetPhi,"fTMetPhi/F");
  iTree.Branch("trimsumet"  ,&fTSumet ,"fTSumet/F");
  iTree.Branch("puppimet"   ,&fPMet    ,"fMet/F");
  iTree.Branch("puppimetphi",&fPMetPhi ,"fMetPhi/F");
  iTree.Branch("puppisumet" ,&fPSumet  ,"fSumet/F");
}
//void puppiJetAnalyzer(
int main( int argc, char **argv ) {
		      std::string iName="/Users/Phil//PUPPI/samples/Zj_80.dat";//) {
    gROOT->ProcessLine("#include <vector>");                
    int nEvts = 0;

    int maxEvents = atoi(argv[1]);

    std::cout << "Processing " << iName << std::endl;
    fin.open(iName.c_str());
    
    TFile* lFile = new TFile("output/OutputMet.root", "RECREATE");
    TTree* lTree      = new TTree("tree", "tree");
    addBranches(*lTree);
 
    std::vector < fastjet::PseudoJet > allParticles;
    std::vector < int > v_isPU;
    std::vector < int > v_isCh;
    readEvent( allParticles, v_isPU, v_isCh );        
    puppiTMVAContainer curEvent(allParticles, v_isPU, v_isCh); //TMVA
   while(true){
      nEvts++;
      if (nEvts % 10 == 0) std::cout << "event no. = " << nEvts << std::endl;
      if (nEvts == maxEvents){ break; }
      if(fin.eof()) break;
      readEvent( allParticles, v_isPU, v_isCh );        
      //puppiCleanContainer curCleanEvent(allParticles, v_isPU, v_isCh);
      curEvent.refresh(allParticles, v_isPU, v_isCh); //TMVA
      std::vector<fastjet::PseudoJet> pfParticles    = curEvent.pfParticles();
      std::vector<fastjet::PseudoJet> genParticles   = curEvent.genParticles();
      std::vector<fastjet::PseudoJet> puppiParticles = curEvent.puppiEvent(7,0.5);
      std::vector<fastjet::PseudoJet> trimParticles;//  = curCleanEvent.trimEvent();
      analyzeEvent(pfParticles,trimParticles,genParticles,puppiParticles );        
      lTree->Fill();
      allParticles.clear();
      v_isPU.clear();
      v_isCh.clear();
    }
    lFile->cd();
    lTree->Write();
}
void readEvent( std::vector< fastjet::PseudoJet > &allParticles, std::vector<int> &v_isPU, std::vector<int> &v_isCh ){
  float npart, px, py, pz, e, pdgid, isCh, isPU = 0;
  int pId = 0;
  while(true){
    if(fin.eof()) break;
    fin  >> npart >> px >> py >> pz >> e >> pdgid >> isCh >> isPU;         
    if (px == 0 && py == 0 && pz == 0 && e == 0) return;
    // fill vector of pseudojets
    fastjet::PseudoJet curPseudoJet( px, py, pz, e );
    if (fabs(curPseudoJet.eta()) < 5){
      int lId = 0; if(isPU) lId++; if(isCh) lId+=2;
      curPseudoJet.set_user_index(lId);
      allParticles.push_back( curPseudoJet );
      v_isPU.push_back(isPU);
      v_isCh.push_back(isCh);            
      pId++;
    }
  }
}
float *met( std::vector< fastjet::PseudoJet > &iParticles) { 
  float *lMet = new float[3]; 
  lMet[0] = 0; 
  lMet[1] = 0; 
  lMet[2] = 0; 
  PseudoJet flow;
  for(unsigned int i=0; i<iParticles.size(); i++){
    flow    += iParticles[i];
    lMet[2] += iParticles[i].pt();
  }
  lMet[0] = flow.pt();
  lMet[1] = flow.phi();
  return lMet;
}
void analyzeEvent(std::vector < fastjet::PseudoJet > &constits, std::vector < fastjet::PseudoJet > &trimconstits, std::vector < fastjet::PseudoJet > &genconstits,std::vector < fastjet::PseudoJet > &puppiconstits) { 
  float *pGenMet     =  met(genconstits);
  float *pPFMet      =  met(constits);
  float *pTrimmedMet =  met(trimconstits);
  float *pPuppiMet   =  met(puppiconstits);
      
  fMet      = pPFMet[0]; 
  fMetPhi   = pPFMet[1]; 
  fSumet    = pPFMet[2]; 
  fGMet     = pGenMet[0]; 
  fGMetPhi  = pGenMet[1]; 
  fGSumet   = pGenMet[2]; 
  fTMet     = pTrimmedMet[0]; 
  fTMetPhi  = pTrimmedMet[1]; 
  fTSumet   = pTrimmedMet[2]; 
  fPMet      = pPuppiMet[0]; 
  fPMetPhi   = pPuppiMet[1]; 
  fPSumet    = pPuppiMet[2]; 
} 

