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
#include "RecoObj.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

ifstream fin;
/////////////////////////////////////////////////////////////////////
// Tree variables
/////////////////////////////////////////////////////////////////////
int   fCount      = 0;
int   fGCount     = 0;

float fPartProb   = 0;
float fPartChi2   = 0; 
float fPartChi2PU = 0;
float fPartPt     = 0; 
float fPartEta    = 0; 
float fPartPhi    = 0; 
float fPartM      = 0;
int   fPartPFType = 0; 
float fPartDepth  = 0; 
float fPartTime   = 0; 
float fPartDZ     = 0; 
float fPartD0     = 0; 

int   fVtxId      = 0; 
float fPartTk     = 0; 
float fPartVtx    = 0; 
float fGXPt       = 0;
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
float fRho    = 0; 

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
void plotEvent( std::vector < fastjet::PseudoJet > constits, std::string name, std::vector < fastjet::PseudoJet > jets );
void setupCMSReadOut(TTree *iTree );
void setupGenCMSReadOut(TTree *iTree );
void readGenCMSEvent(TTree *iTree, std::vector< fastjet::PseudoJet > &allParticles );
void readCMSEvent   (TTree *iTree, std::vector< RecoObj >            &allParticles );
void getJets(std::vector < fastjet::PseudoJet > &constits,std::vector < fastjet::PseudoJet > &jets);
void getCleanJets(std::vector < fastjet::PseudoJet > &constits,std::vector < fastjet::PseudoJet > &jets);
double deltaR(PseudoJet &iJet0,PseudoJet &iJet1);
float *met( std::vector< fastjet::PseudoJet > &iParticles);
void analyzeMETEvent(std::vector < fastjet::PseudoJet > constits, std::vector < fastjet::PseudoJet > trimconstits, std::vector < fastjet::PseudoJet > genconstits,std::vector < fastjet::PseudoJet > puppiconstits);
void setTDRStyle();
void analyzeEvent(TTree *iTree, std::vector < fastjet::PseudoJet > constits,std::vector < fastjet::PseudoJet > chsconstits, std::vector < fastjet::PseudoJet > trimconstits, 
		  std::vector < fastjet::PseudoJet > genconstits, std::vector < fastjet::PseudoJet > gen2constits);
void addBranches( TTree &iTree){
  iTree.Branch("rho",&fRho,"fRho/F");
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

  iTree.Branch("genmatchpt" ,&fClPt  ,"fClPt/F");
  iTree.Branch("genmatcheta",&fClEta ,"fClEta/F");
  iTree.Branch("genmatchphi",&fClPhi ,"fClPhi/F");
  iTree.Branch("genmatchm"  ,&fClM   ,"fClM/F");
  iTree.Branch("genmatchmtr",&fTrClM ,"fTrClM/F");
}
void addMETBranches( TTree &iTree){
  iTree.Branch("rho"      ,&fRho     ,"fRho/F");
  iTree.Branch("genmet"   ,&fGMet    ,"fGMet/F");
  iTree.Branch("genmetphi",&fGMetPhi ,"fGMetPhi/F");
  iTree.Branch("gensumet" ,&fGSumet  ,"fGSumet/F");
  iTree.Branch("met"      ,&fMet     ,"fMet/F");
  iTree.Branch("metphi"   ,&fMetPhi  ,"fMetPhi/F");
  iTree.Branch("sumet"    ,&fSumet   ,"fSumet/F");
  iTree.Branch("cleanmet"  ,&fTMet    ,"fTMet/F");
  iTree.Branch("cleanmetphi",&fTMetPhi,"fTMetPhi/F");
  iTree.Branch("cleansumet" ,&fTSumet ,"fTSumet/F");
  iTree.Branch("puppimet"   ,&fPMet    ,"fMet/F");
  iTree.Branch("puppimetphi",&fPMetPhi ,"fMetPhi/F");
  iTree.Branch("puppisumet" ,&fPSumet  ,"fSumet/F");
}
bool findMuon(std::vector < RecoObj > allParticles ) { 
  for(int i0 = 0; i0 < allParticles.size(); i0++) { 
    if(allParticles[i0].pfType != 3) continue;
    if(allParticles[i0].pt     > 20) return true;
  }
  return false;
}
//void puppiJetAnalyzer(
int main( int argc, char **argv ) {
  setTDRStyle();

  gROOT->ProcessLine("#include <vector>");                
  int nEvts = 0;
  
  int maxEvents  = atoi(argv[1]);
  char *recoFile = argv[2];;
  char *genFile  = argv[3];
  int useMVA     = atoi(argv[4]);
  
  TFile* lReadout = new TFile(recoFile);
  TTree* fTree    = (TTree*) lReadout->FindObjectAny("Result");
  setupCMSReadOut(fTree);
  
  TFile* lGenReadout = new TFile("GenOutput.root");
  TTree* fGenTree    = (TTree*) lGenReadout->FindObjectAny("Result");
  setupGenCMSReadOut(fGenTree);
  
  TFile* lFile = new TFile("output/OutputTmp.root", "RECREATE");
  TTree* lTree      = new TTree("tree", "tree");
  TTree* lMetTree   = new TTree("met", "met");
  addBranches(*lTree);
  addMETBranches(*lMetTree);

  std::vector < RecoObj > allParticles;
  std::vector < fastjet::PseudoJet > genParticles;
  std::vector < fastjet::PseudoJet > pfParticles;
  std::vector < fastjet::PseudoJet > chsParticles;
  std::vector < fastjet::PseudoJet > puppiParticles;
  std::vector < fastjet::PseudoJet > genmatchParticles;
  
  puppiTMVAContainer *mvaEvent = 0; 
  readCMSEvent   (fTree,    allParticles);
  readGenCMSEvent(fGenTree, genParticles);
  if(useMVA) mvaEvent = new puppiTMVAContainer(allParticles);
  if(useMVA) mvaEvent->refresh(allParticles);
  allParticles     .clear();
  genParticles     .clear();
  genmatchParticles.clear();
  pfParticles      .clear();
  chsParticles     .clear();
  puppiParticles   .clear();
  while(true){
      nEvts++;
      if (nEvts % 10 == 0) std::cout << "event no. = " << nEvts << std::endl;
      if (nEvts == maxEvents){ break; }
      readCMSEvent(fTree, allParticles);
      readGenCMSEvent(fGenTree, genParticles);
      if(!useMVA) {
	puppiCleanContainer curEvent(allParticles);
	puppiParticles      = curEvent.puppiEvent(7,0.5);	
        pfParticles         = curEvent.pfParticles();
        chsParticles        = curEvent.pfchsParticles();
        genmatchParticles   = curEvent.genParticles();
      }
      if(useMVA) {
	mvaEvent->refresh(allParticles);
	puppiParticles      = mvaEvent->puppiEvent(7,0.5);	
        pfParticles         = mvaEvent->pfParticles();
        chsParticles        = mvaEvent->pfchsParticles();
        genmatchParticles   = mvaEvent->genParticles();
      }
      analyzeEvent   (lTree,pfParticles,chsParticles,  puppiParticles,genParticles,genmatchParticles);        
      analyzeMETEvent(      pfParticles,puppiParticles,genParticles  ,genmatchParticles);       
      GridMedianBackgroundEstimator lGrid(5.0,0.8);
      lGrid.set_particles(pfParticles);
      fRho = lGrid.rho(); 
      if(nEvts < -10) { 
	std::vector < fastjet::PseudoJet > jets;  
	std::vector < fastjet::PseudoJet > chsjets;  
	std::vector < fastjet::PseudoJet > puppijets;  
	std::vector < fastjet::PseudoJet > genjets;  
	//Get All Jet Colelctions
	getJets(genParticles,    genjets);
	getJets(pfParticles,        jets);
	getJets(chsParticles,    chsjets);
	getJets(puppiParticles,puppijets);
	std::stringstream pSS; pSS << "Label" << nEvts;
	plotEvent(genParticles   ,"gen"  +pSS.str(),genjets );
	plotEvent(pfParticles    ,"pf"   +pSS.str(),jets );
	plotEvent(chsParticles   ,"pfchs"+pSS.str(),chsjets );
	plotEvent(puppiParticles ,"PUPPI"+pSS.str(),puppijets );
      }
      genParticles.clear();
      allParticles.clear();
      pfParticles .clear();
      chsParticles.clear();
      puppiParticles.clear();
      lMetTree->Fill();
    }
    lFile->cd();
    lMetTree->Write();
    lTree   ->Write();
}
void setupCMSReadOut(TTree *iTree ) { 
  fCount = 0;
  iTree->SetBranchAddress("pt"    ,&fPartPt);
  iTree->SetBranchAddress("eta"   ,&fPartEta);
  iTree->SetBranchAddress("phi"   ,&fPartPhi);
  iTree->SetBranchAddress("m"     ,&fPartM);
  iTree->SetBranchAddress("genpt" ,&fGXPt);
  iTree->SetBranchAddress("vtxid" ,&fVtxId);
  iTree->SetBranchAddress("vtx"   ,&fPartVtx);
  iTree->SetBranchAddress("trk"   ,&fPartTk);
  iTree->SetBranchAddress("pftype",&fPartPFType);  
  iTree->SetBranchAddress("depth" ,&fPartDepth);
  iTree->SetBranchAddress("time"  ,&fPartTime);
  iTree->SetBranchAddress("prob"  ,&fPartProb);
  iTree->SetBranchAddress("chi2pv",&fPartChi2);
  iTree->SetBranchAddress("chi2pu",&fPartChi2PU);
  iTree->SetBranchAddress("dZ"   ,&fPartDZ);  
  iTree->SetBranchAddress("d0"    ,&fPartD0);
}
void setupGenCMSReadOut(TTree *iTree ) { 
  fCount = 0;
  iTree->SetBranchAddress("genpt"   ,&fGPt);
  iTree->SetBranchAddress("geneta"  ,&fGEta);
  iTree->SetBranchAddress("genhi"   ,&fGPhi);
  iTree->SetBranchAddress("genm"    ,&fGM);
}
void readCMSEvent(TTree *iTree, std::vector< RecoObj > &allParticles) { 
  float px, py, pz, e, pdgid, isCh, isPU,isPV = 0;
  while(true){
    iTree->GetEntry(fCount);
    fCount++;
    if(fPartPt == -1) break;
    isCh   = (fVtxId > -1);//(fPartPFType == 1 || fPartPFType == 2 || fPartPFType == 3) && (fVtxId > -1 || fabs(fPartDZ) < 0.2) ;
    isPV   = (fVtxId == 0);//  || (fabs(fPartDZ) < 0.2 && isCh);// && fabs(fPartD0) < 0.2));
    isPU   = (fGXPt/fPartPt < 0.2);//*(1-TMath::Min(fGXPt,float(1.))));
    TLorentzVector pVec; pVec.SetPtEtaPhiM(fPartPt,fPartEta,fPartPhi,fPartM);
    int lId = 0; if(!isPV) lId++; if(isCh) lId+=2;
    if(!isPU) lId += 4;
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
    curPseudoJet.expProb = fPartProb;
    curPseudoJet.expChi2PU = fPartChi2PU;
    curPseudoJet.expChi2   = fPartChi2;
    curPseudoJet.d0        = fPartD0;
    curPseudoJet.dZ        = fPartDZ;
    allParticles.push_back( curPseudoJet );
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
    int lId = 0; if(isPU) lId++; if(isCh) lId+=2;
    lId += 4;
    curPseudoJet.set_user_index(lId);
    allParticles.push_back( curPseudoJet );
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
void analyzeEvent(TTree *iTree, std::vector < fastjet::PseudoJet > constits,std::vector < fastjet::PseudoJet > chsconstits, std::vector < fastjet::PseudoJet > puppiconstits, std::vector < fastjet::PseudoJet > genconstits,std::vector < fastjet::PseudoJet > genmatchconstits) { 
  std::vector < fastjet::PseudoJet > jets;  
  std::vector < fastjet::PseudoJet > chsjets;  
  std::vector < fastjet::PseudoJet > puppijets;  
  std::vector < fastjet::PseudoJet > genjets;  
  std::vector < fastjet::PseudoJet > genmatchjets;  
  //Get All Jet Colelctions
  getJets(constits,        jets);
  getJets(chsconstits,     chsjets);
  getJets(puppiconstits,   puppijets);
  getJets(genconstits,     genjets);
  getJets(genmatchconstits,genmatchjets);
  //Trimmer
  fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.3), fastjet::SelectorPtFractionMin(0.05)));
  //Rho 
  GridMedianBackgroundEstimator lGrid(5.0,0.8);
  lGrid.set_particles(constits);
  //Rho on the modified PF candiates
  GridMedianBackgroundEstimator lGridTrim(5.0,0.8);
  lGridTrim.set_particles(puppiconstits);
  GridMedianBackgroundEstimator lGridCHS(5.0,0.8);
  std::vector<PseudoJet> chsCentralConstits; for(unsigned int i0 = 0; i0 < chsconstits.size(); i0++) if(fabs(chsconstits[i0].eta()) < 2.5) chsCentralConstits.push_back(chsconstits[i0]);
  lGridCHS.set_particles(chsCentralConstits);
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
    fRho = lGrid.rho();
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
      if(fabs(chsjets[iId].eta()) > 2.5) pCorrJet -= lGrid   .rho() * pArea;
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
    for(unsigned int i1 = 0; i1 < puppijets.size(); i1++) { 
      if(puppijets[i1].pt() < 5) continue;
      double pDR = deltaR(jets[i0],puppijets[i1]);
      if(pDR > 0.2) continue;
      iId = i1;
      break;
    }
    if(iId > -1) { 
      fTPt  = puppijets[iId].pt();
      fTEta = puppijets[iId].eta();
      fTPhi = puppijets[iId].phi();
      fTM   = puppijets[iId].m();
      //Calculate rho again
      PseudoJet pTCorrJet = puppijets[iId];
      pTCorrJet -= lGridTrim.rho() * pArea;
      fTCPt  = pTCorrJet.pt();
      fTCEta = pTCorrJet.eta();
      fTCPhi = pTCorrJet.phi();
      fTCM   = pTCorrJet.m();
      fastjet::PseudoJet pTrimmedJet = (trimmer)(puppijets[iId]);
      fTrTM  = pTrimmedJet.m();
    }
    iId = -1;
    for(unsigned int i1 = 0; i1 < genmatchjets.size(); i1++) { 
      if(genmatchjets[i1].pt() < 5) continue;
      double pDR = deltaR(jets[i0],genmatchjets[i1]);
      if(pDR > 0.2) continue;
      iId = i1;
      break;
    } 
    if(iId > -1) { 
      fClPt  = genmatchjets[iId].pt();
      fClEta = genmatchjets[iId].eta();
      fClPhi = genmatchjets[iId].phi();
      fClM   = genmatchjets[iId].m();
    }
    iTree->Fill();
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
void analyzeMETEvent(std::vector < fastjet::PseudoJet > constits, std::vector < fastjet::PseudoJet > puppiconstits, std::vector < fastjet::PseudoJet > genconstits,std::vector < fastjet::PseudoJet > genmatchconstits) { 
  float *pGenMet     =  met(genconstits);
  float *pPFMet      =  met(constits);
  float *pTrimmedMet =  met(genmatchconstits);
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


void plotEvent( std::vector < fastjet::PseudoJet > constits, std::string iName, std::vector < fastjet::PseudoJet > jets ) { 
    
    double maxCell = 0;
    TH2F* h2d    = new TH2F( "h2d"   ,";#eta;#phi;pT (GeV)",50, -5,5, 50,0,2*TMath::Pi() );
    TH2F* h2d_PU = new TH2F( "h2d_PU",";#eta;#phi;pT (GeV)",50, -5,5, 50,0,2*TMath::Pi() );
    
    for (unsigned int i = 0; i < constits.size(); i++){
      int curBin = h2d->FindBin(constits[i].eta(),constits[i].phi());
      if(constits[i].pt() > maxCell) maxCell = constits[i].pt();
      if(constits[i].user_index() > 3) h2d   ->SetBinContent( curBin, h2d   ->GetBinContent(curBin) + constits[i].pt() );
      if(constits[i].user_index() < 4) h2d_PU->SetBinContent( curBin, h2d_PU->GetBinContent(curBin) + constits[i].pt() );
    }
    
    TCanvas* can = new TCanvas("can","can",800,600);
    h2d->SetMaximum( 10);//1.1*max(h2d->GetMaximum(),h2d_PU->GetMaximum()) );
    h2d->SetMinimum( 1e-1 );
    h2d_PU->SetMaximum( 10);
    h2d_PU->SetMinimum( 1e-1 );
    h2d   ->Draw("COLZ");
    h2d_PU->SetLineWidth( 1 );
    h2d_PU->Draw("BOX SAME");
    //can->SetLogz();
    
    //draw jets
    for (unsigned j = 0; j < jets.size(); j++){
        TEllipse* cir = new TEllipse(jets[j].eta(),jets[j].phi(),0.7,0.7);
        cir->SetFillStyle(0);
        if (jets[j].pt() > 50 && jets[j].pt() < 200){
            cir->SetLineColor( 7 );
            cir->SetLineWidth( 2 );            
        }
        else if (jets[j].pt() > 200){
            cir->SetLineColor( 4 );
            cir->SetLineWidth( 2 );            
        }
        else{
            cir->SetLineColor( 6 );
        }
        
        cir->Draw("sames");
    }
    can->SaveAs((iName+".png").c_str());
    delete can;
    delete h2d;
    delete h2d_PU;
}
void setTDRStyle() {
    TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
    
    // For the canvas:
    tdrStyle->SetCanvasBorderMode(0);
    tdrStyle->SetCanvasColor(kWhite);
    tdrStyle->SetCanvasDefH(750); //Height of canvas
    tdrStyle->SetCanvasDefW(1050); //Width of canvas
    tdrStyle->SetCanvasDefX(0);   //POsition on screen
    tdrStyle->SetCanvasDefY(0);
    
    // For the Pad:
    tdrStyle->SetPadBorderMode(0);
    // tdrStyle->SetPadBorderSize(Width_t size = 1);
    tdrStyle->SetPadColor(kWhite);
    tdrStyle->SetPadGridX(false);
    tdrStyle->SetPadGridY(false);
    tdrStyle->SetGridColor(0);
    tdrStyle->SetGridStyle(3);
    tdrStyle->SetGridWidth(1);
    
    // For the frame:
    tdrStyle->SetFrameBorderMode(0);
    tdrStyle->SetFrameBorderSize(1);
    tdrStyle->SetFrameFillColor(0);
    tdrStyle->SetFrameFillStyle(0);
    tdrStyle->SetFrameLineColor(1);
    tdrStyle->SetFrameLineStyle(1);
    tdrStyle->SetFrameLineWidth(1);
    
    // For the histo:
    // tdrStyle->SetHistFillColor(1);
    // tdrStyle->SetHistFillStyle(0);
    tdrStyle->SetHistLineColor(1);
    tdrStyle->SetHistLineStyle(0);
    tdrStyle->SetHistLineWidth(1);
    // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
    // tdrStyle->SetNumberContours(Int_t number = 20);
    
    tdrStyle->SetEndErrorSize(2);
    //  tdrStyle->SetErrorMarker(20);
    tdrStyle->SetErrorX(0.);
    
    tdrStyle->SetMarkerStyle(20);
    
    //For the fit/function:
    tdrStyle->SetOptFit(1);
    tdrStyle->SetFitFormat("5.4g");
    tdrStyle->SetFuncColor(2);
    tdrStyle->SetFuncStyle(1);
    tdrStyle->SetFuncWidth(1);
    
    //For the date:
    tdrStyle->SetOptDate(0);
    // tdrStyle->SetDateX(Float_t x = 0.01);
    // tdrStyle->SetDateY(Float_t y = 0.01);
    
    // For the statistics box:
    tdrStyle->SetOptFile(0);
    tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
    tdrStyle->SetStatColor(kWhite);
    tdrStyle->SetStatFont(42);
    tdrStyle->SetStatFontSize(0.010);
    tdrStyle->SetStatTextColor(1);
    tdrStyle->SetStatFormat("6.4g");
    tdrStyle->SetStatBorderSize(1);
    tdrStyle->SetStatH(0.25);
    tdrStyle->SetStatW(0.15);
    // tdrStyle->SetStatStyle(Style_t style = 1001);
    // tdrStyle->SetStatX(Float_t x = 0);
    // tdrStyle->SetStatY(Float_t y = 0);
    
    // Margins:
    tdrStyle->SetPadTopMargin(0.05);
    tdrStyle->SetPadBottomMargin(0.13);
    tdrStyle->SetPadLeftMargin(0.17);
    tdrStyle->SetPadRightMargin(0.17);
    
    // For the Global title:
    
    tdrStyle->SetOptTitle(0);
    tdrStyle->SetTitleFont(42);
    tdrStyle->SetTitleColor(1);
    tdrStyle->SetTitleTextColor(1);
    tdrStyle->SetTitleFillColor(10);
    tdrStyle->SetTitleFontSize(0.005);
    // tdrStyle->SetTitleH(0); // Set the height of the title box
    // tdrStyle->SetTitleW(0); // Set the width of the title box
    // tdrStyle->SetTitleX(0); // Set the position of the title box
    // tdrStyle->SetTitleY(0.985); // Set the position of the title box
    // tdrStyle->SetTitleStyle(Style_t style = 1001);
    // tdrStyle->SetTitleBorderSize(2);
    
    // For the axis titles:
    
    tdrStyle->SetTitleColor(1, "XYZ");
    tdrStyle->SetTitleFont(42, "XYZ");
    tdrStyle->SetTitleSize(0.06, "XYZ");
    // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
    // tdrStyle->SetTitleYSize(Float_t size = 0.02);
    tdrStyle->SetTitleXOffset(0.9);
    tdrStyle->SetTitleYOffset(1.25);
    // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
    
    // For the axis labels:
    
    tdrStyle->SetLabelColor(1, "XYZ");
    tdrStyle->SetLabelFont(42, "XYZ");
    tdrStyle->SetLabelOffset(0.007, "XYZ");
    tdrStyle->SetLabelSize(0.05, "XYZ");
    
    // For the axis:
    
    tdrStyle->SetAxisColor(1, "XYZ");
    tdrStyle->SetStripDecimals(kTRUE);
    tdrStyle->SetTickLength(0.03, "XYZ");
    tdrStyle->SetNdivisions(505, "XYZ");
    tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    tdrStyle->SetPadTickY(1);
    
    // Change for log plots:
    tdrStyle->SetOptLogx(0);
    tdrStyle->SetOptLogy(0);
    tdrStyle->SetOptLogz(0);
    
    // Postscript options:
    tdrStyle->SetPaperSize(20.,20.);
    // tdrStyle->SetLineScalePS(Float_t scale = 3);
    // tdrStyle->SetLineStyleString(Int_t i, const char* text);
    // tdrStyle->SetHeaderPS(const char* header);
    // tdrStyle->SetTitlePS(const char* pstitle);
    
    // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
    // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
    // tdrStyle->SetPaintTextFormat(const char* format = "g");
    // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    // tdrStyle->SetTimeOffset(Double_t toffset);
    // tdrStyle->SetHistMinimumZero(kTRUE);
    
    tdrStyle->SetPalette(1);
    
    
    tdrStyle->cd();
    
}
