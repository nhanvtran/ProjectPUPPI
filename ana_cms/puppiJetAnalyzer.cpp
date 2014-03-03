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

#include <fastjet/tools/GridMedianBackgroundEstimator.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include <fastjet/contrib/JetCleanser.hh>

#include "puppiCleanContainer.hh"
#include "puppiTMVAContainer.hh"
#include "NoTrees.hh"


using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

ifstream fin;
/////////////////////////////////////////////////////////////////////
// Tree variables
/////////////////////////////////////////////////////////////////////
int   fCount      = 0;
int   fGCount     = 0;
float fPartPt     = 0; 
float fPartEta    = 0; 
float fPartPhi    = 0; 
float fPartM      = 0;
float fGXPt       = 0;
float fPartTk     = 0; 
float fPartVtx    = 0; 
int   fVtxId      = 0; 
float fGPt        = 0; 
float fGEta       = 0; 
float fGPhi       = 0; 
float fGM         = 0;

float fPt     = 0; 
float fEta    = 0; 
float fPhi    = 0; 
float fM      = 0;
float fTrM    = 0;
float fCPt    = 0; 
float fCEta   = 0; 
float fCPhi   = 0; 
float fCM     = 0;
float fTrCM   = 0;
float fTPt    = 0; 
float fTEta   = 0; 
float fTPhi   = 0; 
float fTM     = 0; 
float fTrTM   = 0; 
float fTCPt   = 0; 
float fTCEta  = 0; 
float fTCPhi  = 0; 
float fTCM    = 0; 
float fGenPt  = 0; 
float fGenEta = 0; 
float fGenPhi = 0; 
float fGenM   = 0; 
float fTrGenM = 0; 
float fClPt   = 0; 
float fClEta  = 0; 
float fClPhi  = 0; 
float fClM    = 0; 
float fTrClM  = 0;
/////////////////////////////////////////////////////////////////////
// Helper functions
/////////////////////////////////////////////////////////////////////
void setupCMSReadOut   (TTree *iTree );
void setupGenCMSReadOut(TTree *iTree );
void readGenCMSEvent   (TTree *iTree, std::vector< fastjet::PseudoJet > &allParticles );
void readCMSEvent      (TTree *iTree, std::vector< fastjet::PseudoJet > &allParticles, 
			std::vector<int> &v_isPU, std::vector<int> &v_isCh,std::vector<float> &v_pTk,std::vector<float> &v_pVtx,std::vector<int> &v_VId );
void readEvent( std::vector< fastjet::PseudoJet > &allParticles, std::vector<int> &v_isPU, std::vector<int> &v_isCh );
double deltaR(PseudoJet &iJet0,PseudoJet &iJet1);
void getJets(std::vector < fastjet::PseudoJet > &constits,std::vector < fastjet::PseudoJet > &jets);
void getCleanJets(std::vector < fastjet::PseudoJet > &constits,std::vector < fastjet::PseudoJet > &jets);
void analyzeEvent(TTree *iTree, std::vector < fastjet::PseudoJet > constits,std::vector < fastjet::PseudoJet > chsconstits, std::vector < fastjet::PseudoJet > trimconstits, 
		  std::vector < fastjet::PseudoJet > genconstits, std::vector < fastjet::PseudoJet > gen2constits);
void addBranches( TTree &iTree){
  iTree.Branch("pt" ,&fPt ,"fPt/F");
  iTree.Branch("eta",&fEta,"fEta/F");
  iTree.Branch("phi",&fPhi,"fPhi/F");
  iTree.Branch("m"  ,&fM  ,"fM/F");
  iTree.Branch("mtr",&fTrM,"fTrM/F");

  iTree.Branch("chspt" ,&fCPt ,"fCPt/F");
  iTree.Branch("chseta",&fCEta,"fCEta/F");
  iTree.Branch("chsphi",&fCPhi,"fCPhi/F");
  iTree.Branch("chsm"  ,&fCM  ,"fCM/F");
  iTree.Branch("chsmtr",&fTrCM,"fTrCM/F");

  iTree.Branch("trimpt" ,&fTPt ,"fTPt/F");
  iTree.Branch("trimeta",&fTEta,"fTEta/F");
  iTree.Branch("trimphi",&fTPhi,"fTPhi/F");
  iTree.Branch("trimm"  ,&fTM  ,"fTM/F");
  iTree.Branch("trimmtr",&fTrTM,"fTTrM/F");

  iTree.Branch("trimcorrpt" ,&fTCPt ,"fTCPt/F");
  iTree.Branch("trimcorreta",&fTCEta,"fTCEta/F");
  iTree.Branch("trimcorrphi",&fTCPhi,"fTCPhi/F");
  iTree.Branch("trimcorrm"  ,&fTCM  ,"fTCM/F");

  iTree.Branch("genpt" ,&fGenPt ,"fGenPt/F");
  iTree.Branch("geneta",&fGenEta,"fGenEta/F");
  iTree.Branch("genphi",&fGenPhi,"fGenPhi/F");
  iTree.Branch("genm"  ,&fGenM  ,"fGenM/F");
  iTree.Branch("genmtr",&fTrGenM ,"fTrGenM/F");

  iTree.Branch("cleanpt" ,&fClPt  ,"fClPt/F");
  iTree.Branch("cleaneta",&fClEta ,"fClEta/F");
  iTree.Branch("cleanphi",&fClPhi ,"fClPhi/F");
  iTree.Branch("cleanm"  ,&fClM   ,"fClM/F");
  iTree.Branch("cleanmtr",&fTrClM ,"fTrClM/F");
}
//void puppiJetAnalyzer(
int main( int argc, char **argv ) {
  std::string iName="/Users/Phil//PUPPI/samples/Zj_80.dat";//) {
  gROOT->ProcessLine("#include <vector>");                
  int nEvts = 0;
  
  int maxEvents = atoi(argv[1]);
  
  std::cout << "Processing " << iName << std::endl;
  TFile* lReadout = new TFile("RecoOutput.root");
  TTree* fTree    = (TTree*) lReadout->FindObjectAny("Result");
  setupCMSReadOut(fTree);
  
  TFile* lGenReadout = new TFile("GenOutput.root");
  TTree* fGenTree    = (TTree*) lGenReadout->FindObjectAny("Result");
  setupGenCMSReadOut(fGenTree);
  
  TFile* lFile = new TFile("output/OutputTmp.root", "RECREATE");
  TTree* lTree      = new TTree("tree", "tree");
  addBranches(*lTree);
  
  std::vector < fastjet::PseudoJet > allParticles;
  std::vector < fastjet::PseudoJet > genParticles;
  std::vector < int > v_isPU;
  std::vector < int > v_isCh;
  std::vector < float > v_Tk;
  std::vector < float > v_Vtx;
  std::vector < int >   v_VId;
  
  //puppiTMVAContainer curEvent(allParticles, v_isPU, v_isCh,v_Tk,v_Vtx,v_VId); //TMVA
  readCMSEvent   (fTree,    allParticles, v_isPU, v_isCh,v_Tk,v_Vtx,v_VId);        
  readGenCMSEvent(fGenTree, genParticles);
  //curEvent.refresh(allParticles, v_isPU, v_isCh,v_Tk,v_Vtx,v_VId); //TMVA
  genParticles.clear();
  allParticles.clear();
  while(true){
    nEvts++;
    if (nEvts % 10 == 0) std::cout << "event no. = " << nEvts << std::endl;
    if (nEvts == maxEvents){ break; }
    if(fin.eof()) break;
    readCMSEvent(fTree, allParticles, v_isPU, v_isCh,v_Tk,v_Vtx,v_VId);        
    readGenCMSEvent(fGenTree, genParticles);
    puppiCleanContainer curEvent(allParticles, v_isPU, v_isCh);
    //curEvent.refresh(allParticles, v_isPU, v_isCh,v_Tk,v_Vtx,v_VId); //TMVA
    std::vector<fastjet::PseudoJet> puppiParticles = curEvent.puppiEvent(7,0.5);
    analyzeEvent(lTree,curEvent.pfParticles(),curEvent.pfchsParticles(),puppiParticles,genParticles,curEvent.genParticles());        
    genParticles.clear();
    allParticles.clear();
    v_isPU.clear();
    v_isCh.clear();
    v_Tk.clear();
    v_Vtx.clear();
    v_VId.clear();
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
  iTree->SetBranchAddress("genpt",&fGXPt);
  iTree->SetBranchAddress("vtxid",&fVtxId);
  iTree->SetBranchAddress("vtx"  ,&fPartVtx);
  iTree->SetBranchAddress("trk"  ,&fPartTk);
}
void setupGenCMSReadOut(TTree *iTree ) { 
  fCount = 0;
  iTree->SetBranchAddress("genpt"   ,&fGPt);
  iTree->SetBranchAddress("geneta"  ,&fGEta);
  iTree->SetBranchAddress("genhi"   ,&fGPhi);
  iTree->SetBranchAddress("genm"    ,&fGM);
}
void readCMSEvent(TTree *iTree, std::vector< fastjet::PseudoJet > &allParticles, std::vector<int> &v_isPU, std::vector<int> &v_isCh,
		  std::vector<float> &v_pTk,std::vector<float> &v_pVtx,std::vector<int> &v_pVId ){
  float px, py, pz, e, pdgid, isCh, isPU, isPUGen = 0;
  while(true){
    iTree->GetEntry(fCount);
    fCount++;
    if(fPartPt == -1) break;
    isCh = (fVtxId != -1);
    isPUGen = (fGXPt/fPartPt < 0.3);
    isPU = (fVtxId > 0);
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
      if(isPUGen) lId += 4;
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
void readGenCMSEvent(TTree *iTree, std::vector< fastjet::PseudoJet > &allParticles) { 
  float px, py, pz, e, pdgid, isCh, isPU = 0;
  while(true){
    iTree->GetEntry(fGCount);
    fGCount++;
    if(fGPt == -1) break;
    TLorentzVector pVec; pVec.SetPtEtaPhiM(fGPt,fGEta,fGPhi,fGM);
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
    }
  }
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
      int lId = 0; if(isPU) lId++; 
      if(isCh && fabs(curPseudoJet.eta()) < 2.5) lId+=2;
      curPseudoJet.set_user_index(lId);
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
  return;
}
void getCleanJets(std::vector < fastjet::PseudoJet > &constits,std::vector < fastjet::PseudoJet > &jets) { 
  double rParam = 0.7;
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);    
  // do cleansing
  bool doCleansing = false; 
  JetCleanser linear_cleanser_B(0.25, JetCleanser::linear_cleansing, JetCleanser::input_nc_separate);
  linear_cleanser_B.SetLinearParameters(0.65);
  vector<PseudoJet> p_chLV;
  vector<PseudoJet> p_chPU;
  vector<PseudoJet> p_neut;    
  for (unsigned int i = 0; i < constits.size(); i++){
    if      (constits[i].user_index() == 1) p_neut.push_back(constits[i]);
    else if (constits[i].user_index() == 2) p_chLV.push_back(constits[i]);
    else if (constits[i].user_index() == 3) p_chPU.push_back(constits[i]);        
    else continue;
  }
  if (p_chPU.size() > 0) doCleansing = true;
  
  if (doCleansing){
    vector< vector<fastjet::PseudoJet> > sets;
    sets.push_back( constits );           // calorimeter cells
    sets.push_back( p_chLV );             // tracks from primary interaction
    sets.push_back( p_chPU );             // tracks from pileup
    sets.push_back( p_neut );             // neutral particles
    
    // collect jets
    vector< vector<fastjet::PseudoJet> > jet_sets = ClusterSets(jetDef, constits, sets, 25.0);
    vector<fastjet::PseudoJet> jets_plain     = jet_sets[0];
    vector<fastjet::PseudoJet> jets_tracks_LV = jet_sets[1];
    vector<fastjet::PseudoJet> jets_tracks_PU = jet_sets[2];
    vector<fastjet::PseudoJet> jets_neutrals  = jet_sets[3];
    
    for (unsigned int i=0; i<jets_plain.size(); i++){
      PseudoJet plain_jet = jets_plain[i];
      PseudoJet lin_cleansed_jet = linear_cleanser_B( jets_neutrals[i].constituents(), jets_tracks_LV[i].constituents(), jets_tracks_PU[i].constituents() );
      jets.push_back(lin_cleansed_jet);
    }
  }   
}
double deltaR(PseudoJet &iJet0,PseudoJet &iJet1) { 
  double pDPhi = iJet0.phi()-iJet1.phi();
  double pDEta = iJet0.eta()-iJet1.eta();
  if(fabs(pDPhi) > 2.*TMath::Pi()-fabs(pDPhi)) pDPhi =  2.*TMath::Pi()-fabs(pDPhi);
  return sqrt(pDPhi*pDPhi+pDEta*pDEta);
}
void analyzeEvent(TTree *iTree, std::vector < fastjet::PseudoJet > constits,std::vector < fastjet::PseudoJet > chsconstits, std::vector < fastjet::PseudoJet > trimconstits, std::vector < fastjet::PseudoJet > genconstits,std::vector < fastjet::PseudoJet > gen2constits) { 
  std::vector < fastjet::PseudoJet > jets;  
  std::vector < fastjet::PseudoJet > cleanjets;  
  std::vector < fastjet::PseudoJet > trimjets;  
  std::vector < fastjet::PseudoJet > genjets;  
  std::vector < fastjet::PseudoJet > chsjets;  
  //Get All Jet Colelctions
  getJets(constits,jets);
  getJets(gen2constits,cleanjets);
  getJets(trimconstits,trimjets);
  getJets(genconstits,genjets);
  getJets(chsconstits,chsjets);
  //Trimmer
  fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.3), fastjet::SelectorPtFractionMin(0.05)));
  //Rho 
  GridMedianBackgroundEstimator lGrid(5.0,0.8);
  lGrid.set_particles(constits);
  //Rho on the modified PF candiates
  GridMedianBackgroundEstimator lGridTrim(5.0,0.8);
  lGridTrim.set_particles(trimconstits);
  GridMedianBackgroundEstimator lGridCHS(5.0,0.8);
  lGridCHS.set_particles(chsconstits);
  for(unsigned int i0 = 0; i0 < jets.size(); i0++) { 
    fPt    = -20; fEta    = -20; fPhi    = -20; fM    = -20; fTrM    = -20;
    fCPt   = -20; fCEta   = -20; fCPhi   = -20; fCM   = -20; fTrCM   = -20; 
    fGenPt = -20; fGenEta = -20; fGenPhi = -20; fGenM = -20; fTrGenM = -20;
    fTPt   = -20; fTEta   = -20; fTPhi   = -20; fTM   = -20; fTrTM   = -20;
    fTCPt  = -20; fTCEta  = -20; fTCPhi  = -20; fTCM  = -20; 
    fClPt  = -20; fClEta  = -20; fClPhi  = -20; fClM  = -20; fTrClM  = -20;
    fastjet::PseudoJet pTTrimmedJet = (trimmer)(jets[i0]);
    //Corrected PF Jets
    PseudoJet pCorrJet = jets[i0];
    PseudoJet pArea    = jets[i0].area_4vector();
    //double    pTArea   = pTTrimmedJet.area();
    pCorrJet     -= lGrid.rho() * pArea;
    fPt  = pCorrJet.pt();
    fEta = pCorrJet.eta();
    fPhi = pCorrJet.phi();
    fM   = pCorrJet.m();
    fTrM = pTTrimmedJet.m();
    //Match to CHS Jets
    int iId = -1;
    for(unsigned int i1 = 0; i1 < chsjets.size(); i1++) { 
      if(chsjets[i1].pt() < 5) continue;
      double pDR = deltaR(jets[i0],chsjets[i1]);
      if(pDR > 0.2) continue;
      iId = i1;
      break;
    }
    if(iId > -1) { 
      pCorrJet  = chsjets[iId];
      pArea     = chsjets[iId].area_4vector();
      pCorrJet -= lGridCHS.rho() * pArea;
      pTTrimmedJet = (trimmer)(chsjets[iId]);
      fCPt  = pCorrJet.pt();
      fCEta = pCorrJet.eta();
      fCPhi = pCorrJet.phi();
      fCM   = pCorrJet.m();
      fTrCM = pTTrimmedJet.m();
    }
    iId = -1;
    //Match to Gen Jets
    for(unsigned int i1 = 0; i1 < genjets.size(); i1++) { 
      if(genjets[i1].pt() < 5) continue;
      double pDR = deltaR(jets[i0],genjets[i1]);
      if(pDR > 0.2) continue;
      iId = i1;
      break;
    }
    if(iId > -1) { 
      fGenPt  = genjets[iId].pt();
      fGenEta = genjets[iId].eta();
      fGenPhi = genjets[iId].phi();
      fGenM   = genjets[iId].m();
      fastjet::PseudoJet pTrimmedJet = (trimmer)(genjets[iId]);
      fTrGenM = pTrimmedJet.m();
    }
    iId = -1;
    for(unsigned int i1 = 0; i1 < trimjets.size(); i1++) { 
      if(trimjets[i1].pt() < 5) continue;
      double pDR = deltaR(jets[i0],trimjets[i1]);
      if(pDR > 0.2) continue;
      iId = i1;
      break;
    }
    if(iId > -1) { 
      fTPt  = trimjets[iId].pt();
      fTEta = trimjets[iId].eta();
      fTPhi = trimjets[iId].phi();
      fTM   = trimjets[iId].m();
      //Calculate rho again
      PseudoJet pTCorrJet = trimjets[iId];
      //PseudoJet pTArea    = trimjets[iId].area_4vector();
      pTCorrJet -= lGridTrim.rho() * pArea;
      fTCPt  = pTCorrJet.pt();
      fTCEta = pTCorrJet.eta();
      fTCPhi = pTCorrJet.phi();
      fTCM   = pTCorrJet.m();
      fastjet::PseudoJet pTrimmedJet = (trimmer)(trimjets[iId]);
      fTrTM  = pTrimmedJet.m();
    }
    iId = -1;
    for(unsigned int i1 = 0; i1 < cleanjets.size(); i1++) { 
      if(cleanjets[i1].pt() < 5) continue;
      double pDR = deltaR(jets[i0],cleanjets[i1]);
      if(pDR > 0.2) continue;
      iId = i1;
      break;
    } 
    if(iId > -1) { 
      fClPt  = cleanjets[iId].pt();
      fClEta = cleanjets[iId].eta();
      fClPhi = cleanjets[iId].phi();
      fClM   = cleanjets[iId].m();
      //fastjet::PseudoJet pTrimmedJet = (trimmer)(cleanjets[iId]);
      //fTrClM  = pTrimmedJet.m();
    }
    iTree->Fill();
  }
}
