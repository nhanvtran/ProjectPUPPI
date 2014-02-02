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

#include <fastjet/GridMedianBackgroundEstimator.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

#include "puppiCleanContainer.hh"
#include "NoTrees.hh"

using namespace std;
using namespace fastjet;

ifstream fin;
/////////////////////////////////////////////////////////////////////
// Tree variables
/////////////////////////////////////////////////////////////////////
float fPt     = 0; 
float fEta    = 0; 
float fPhi    = 0; 
float fM      = 0;
float fCPt    = 0; 
float fCEta   = 0; 
float fCPhi   = 0; 
float fCM     = 0;
float fTPt    = 0; 
float fTEta   = 0; 
float fTPhi   = 0; 
float fTM     = 0; 
float fGenPt  = 0; 
float fGenEta = 0; 
float fGenPhi = 0; 
float fGenM   = 0; 
/////////////////////////////////////////////////////////////////////
// Helper functions
/////////////////////////////////////////////////////////////////////
void readEvent( std::vector< fastjet::PseudoJet > &allParticles, std::vector<int> &v_isPU, std::vector<int> &v_isCh );
void getJets(std::vector < fastjet::PseudoJet > &constits,std::vector < fastjet::PseudoJet > &jets);
void analyzeEvent(TTree *iTree, std::vector < fastjet::PseudoJet > constits, std::vector < fastjet::PseudoJet > trimconstits, std::vector < fastjet::PseudoJet > genconstits);
void addBranches( TTree &iTree){
  iTree.Branch("pt" ,&fPt ,"fPt/F");
  iTree.Branch("eta",&fEta,"fEta/F");
  iTree.Branch("phi",&fPhi,"fPhi/F");
  iTree.Branch("m"  ,&fM  ,"fM/F");

  iTree.Branch("corrpt" ,&fCPt ,"fCPt/F");
  iTree.Branch("correta",&fCEta,"fCEta/F");
  iTree.Branch("corrphi",&fCPhi,"fCPhi/F");
  iTree.Branch("corrm"  ,&fCM  ,"fCM/F");

  iTree.Branch("trimpt" ,&fTPt ,"fTPt/F");
  iTree.Branch("trimeta",&fTEta,"fTEta/F");
  iTree.Branch("trimphi",&fTPhi,"fTPhi/F");
  iTree.Branch("trimm"  ,&fTM  ,"fTM/F");

  iTree.Branch("genpt" ,&fGenPt ,"fGenPt/F");
  iTree.Branch("geneta",&fGenEta,"fGenEta/F");
  iTree.Branch("genphi",&fGenPhi,"fGenPhi/F");
  iTree.Branch("genm"  ,&fGenM  ,"fGenM/F");
}
void puppiJetAnalyzer(std::string iName="/Users/Phil//PUPPI/samples/Zj_80.dat") {
    gROOT->ProcessLine("#include <vector>");                
    int nEvts = 0;
    int maxEvents = 100;

    std::cout << "Processing " << iName << std::endl;
    fin.open(iName.c_str());
    
    TFile* lFile = new TFile("output/OutputTmp.root", "RECREATE");
    TTree* lTree      = new TTree("tree", "tree");
    addBranches(*lTree);
 
    std::vector < fastjet::PseudoJet > allParticles;
    std::vector < int > v_isPU;
    std::vector < int > v_isCh;
    readEvent( allParticles, v_isPU, v_isCh );        
    while(true){
      nEvts++;
      if (nEvts % 10 == 0) std::cout << "event no. = " << nEvts << std::endl;
      if (nEvts == maxEvents){ break; }
      if(fin.eof()) break;
      readEvent( allParticles, v_isPU, v_isCh );        
      puppiCleanContainer curEvent(allParticles, v_isPU, v_isCh);
      std::vector<fastjet::PseudoJet> puppiParticles = curEvent.puppiEvent(7,0.5);
      analyzeEvent(lTree,curEvent.pfParticles(),puppiParticles,curEvent.genParticles() );        
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
      curPseudoJet.set_user_index(pId);
      allParticles.push_back( curPseudoJet );
      v_isPU.push_back(isPU);
      v_isCh.push_back(isCh);            
      pId++;
    }
  }
}
void getJets(std::vector < fastjet::PseudoJet > &constits,std::vector < fastjet::PseudoJet > &jets) { 

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
  //double lArea = thisClustering_->area(jets[0]); 
  //fastjet::PseudoJet lAreaVec = thisClustering_->area_4vector(jets[0]); 
  //cout << "Area : " << lArea << " -- " << lAreaVec.pt() << endl;
  return;
  //Selector rapmax = SelectorAbsRapMax(5.0);
  //Subtractor                    lSub (&lGrid);
  //double lRho                  = thisClustering_->median_pt_per_unit_area(rapmax); 
  //double lRho4V                = thisClustering_->median_pt_per_unit_area_4vector(rapmax);
  //double lPtCorr               = out_jets[0].pt() - lRho * lArea; // or:
  //fastjet::PseudoJet lPtCorr4V = out_jets[0] - lRho4V * lAreaVec; 
  //cout << "rho " << lGrid.rho() << " -- " << lRho << " -- rho 4v:  " << lRho4V << " -- pt :  " << lPtCorr << " -- pt corr 4v :  " << lPtCorr4V.pt() << endl;
}
void analyzeEvent(TTree *iTree, std::vector < fastjet::PseudoJet > constits, std::vector < fastjet::PseudoJet > trimconstits, std::vector < fastjet::PseudoJet > genconstits) { 
  std::vector < fastjet::PseudoJet > jets;  
  std::vector < fastjet::PseudoJet > trimjets;  
  std::vector < fastjet::PseudoJet > genjets;  
  //Get All Jet Colelctions
  getJets(constits,jets);
  getJets(trimconstits,trimjets);
  getJets(genconstits,genjets);
  //Rho 
  GridMedianBackgroundEstimator lGrid(5.0,0.8);
  lGrid.set_particles(constits);
  for(unsigned int i0 = 0; i0 < jets.size(); i0++) { 
    fPt    = -20; fEta    = -20; fPhi    = -20; fM    = -20; 
    fCPt   = -20; fCEta   = -20; fCPhi   = -20; fCM   = -20; 
    fGenPt = -20; fGenEta = -20; fGenPhi = -20; fGenM = -20; 
    fTPt   = -20; fTEta   = -20; fTPhi   = -20; fTM   = -20; 
    //Basics
    fPt  = jets[i0].pt();
    fEta = jets[i0].eta();
    fPhi = jets[i0].phi();
    fM   = jets[i0].m();
    //Corrected PF Jets
    PseudoJet pCorrJet = jets[i0];
    PseudoJet pArea    = jets[i0].area_4vector();
    pCorrJet -= lGrid.rho() * pArea;
    fCPt  = pCorrJet.pt();
    fCEta = pCorrJet.eta();
    fCPhi = pCorrJet.phi();
    fCM   = pCorrJet.m();
    int iId = -1;
    //Match to Gen Jets
    for(unsigned int i1 = 0; i1 < genjets.size(); i1++) { 
      if(genjets[i0].pt() < 5) continue;
      double pDR = jets[i0].delta_R(genjets[i1]);
      if(pDR > 0.4) continue;
      iId = i1;
      break;
    }
    if(iId > -1) { 
      fGenPt  = genjets[iId].pt();
      fGenEta = genjets[iId].eta();
      fGenPhi = genjets[iId].phi();
      fGenM   = genjets[iId].m();
    }
    iId = -1;
    for(unsigned int i1 = 0; i1 < trimjets.size(); i1++) { 
      if(trimjets[i0].pt() < 5) continue;
      double pDR = jets[i0].delta_R(trimjets[i1]);
      if(pDR > 0.4) continue;
      iId = i1;
      break;
    }
    if(iId > -1) { 
      fTPt  = trimjets[i0].pt();
      fTEta = trimjets[i0].eta();
      fTPhi = trimjets[i0].phi();
      fTM   = trimjets[i0].m();
    }
    iTree->Fill();
  }
}
