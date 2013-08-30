#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
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

#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

#include "NoTrees.hh"

using namespace std;
using namespace fastjet;


////////////////////-----------------------------------------------
////////////////////-----------------------------------------------
// Global variables

// read in the file 
ifstream fin;

// mass quantities
float jet_m_[10];
float jet_pt_[10];

// Groomed quantities
//     trimmed1 = 
//     pruned1  = 
float jet_m_trimmed1_[10];
float jet_m_pruned1_[10];
float jet_pt_trimmed1_[10];
float jet_pt_pruned1_[10];
float jet_area_[10];
float jet_area_trimmed1_[10];
float jet_area_pruned1_[10];
float jet_nconstituents_[10];

//event
float jet_njets_;

float jet_HT_, jet_PTMiss_;
float jet_HT_trimmed1_, jet_PTMiss_trimmed1_;
float prt_HT_;
float prt_PTMiss_;

float jwj_Nj_;
float jwj_Nj_trimmed1_;
float jwj_HT_;
float jwj_HT_trimmed1_;
float jwj_PTMiss_;
float jwj_PTMiss_trimmed1_;


///////////////////per particle
float pt_part;
float eta_part;

void readEvent( std::vector< fastjet::PseudoJet > &gen, std::vector< fastjet::PseudoJet > &pfchs, std::vector< fastjet::PseudoJet > &pf );
//void clusterEvent( std::vector < fastjet::PseudoJet > constits, std::vector< float > pdgids, TTree &tree );
void analyzeEvent( std::vector < fastjet::PseudoJet > constits, TTree &tree );
void plotEvent( std::vector < fastjet::PseudoJet > genParticles, char* name );
void setTDRStyle();

void initVars(){
    for (int i = 0; i < 10; i++ ){
        jet_m_[i] = -1.;
        jet_m_trimmed1_[i] = -1;
        jet_m_pruned1_[i] = -1;        
        jet_pt_[i] = -1.;   
        jet_pt_trimmed1_[i] = -1;
        jet_pt_pruned1_[i] = -1;        
        jet_area_[i] = -1.;   
        jet_area_trimmed1_[i] = -1;
        jet_area_pruned1_[i] = -1;        
        jet_nconstituents_[i] = -1.;
    }
    jet_njets_ = -1.;
    
    jet_HT_ = -1.;
    jet_PTMiss_ = -1.;
    jet_HT_trimmed1_ = -1.;
    jet_PTMiss_trimmed1_ = -1.;
    prt_HT_ = -1.;
    prt_PTMiss_ = -1.;
    
    jwj_Nj_ = -1.;
    jwj_Nj_trimmed1_ = -1.;
    jwj_HT_ = -1.;
    jwj_HT_trimmed1_ = -1.;
    jwj_PTMiss_ = -1.;
    jwj_PTMiss_trimmed1_ = -1.;

}

void addBranches( TTree &tree ){
    
    tree.Branch("jet_m",&jet_m_,"mjet_[10]/F");
    tree.Branch("jet_m_trimmed1",&jet_m_trimmed1_,"jet_m_trimmed1[10]/F");
    tree.Branch("jet_m_pruned1",&jet_m_pruned1_,"jet_m_pruned1[10]/F");    
    
    tree.Branch("jet_pt",&jet_pt_,"pt[10]/F");
    tree.Branch("jet_pt_trimmed1",&jet_pt_trimmed1_,"jet_pt_trimmed1[10]/F");
    tree.Branch("jet_pt_pruned1",&jet_pt_pruned1_,"jet_pt_pruned1[10]/F");        
    
    tree.Branch("jet_area",&jet_area_,"area[10]/F");
    tree.Branch("jet_area_trimmed1",&jet_area_trimmed1_,"jet_area_trimmed1[10]/F");
    tree.Branch("jet_area_pruned1",&jet_area_pruned1_,"jet_area_pruned1[10]/F");        
    
    tree.Branch("jet_nconstituents",&jet_nconstituents_,"jet_nconstituents[10]/F");
    tree.Branch("jet_njets",&jet_njets_,"jet_njets/F");
        
    tree.Branch("jet_HT",&jet_HT_,"jet_HT/F");
    tree.Branch("jet_PTMiss",&jet_PTMiss_,"jet_PTMiss/F");
    tree.Branch("jet_HT_trimmed1",&jet_HT_trimmed1_,"jet_HT_trimmed1/F");
    tree.Branch("jet_PTMiss_trimmed1",&jet_PTMiss_trimmed1_,"jet_PTMiss_trimmed1/F");
    tree.Branch("prt_HT",&prt_HT_,"prt_HT/F");
    tree.Branch("prt_PTMiss",&prt_PTMiss_,"prt_PTMiss/F");

    tree.Branch("jwj_Nj",&jwj_Nj_,"jwj_Nj/F");
    tree.Branch("jwj_Nj_trimmed1",&jwj_Nj_trimmed1_,"jwj_Nj_trimmed1/F");
    tree.Branch("jwj_HT",&jwj_HT_,"jwj_HT/F");
    tree.Branch("jwj_HT_trimmed1",&jwj_HT_trimmed1_,"jwj_HT_trimmed1/F");
    tree.Branch("jwj_PTMiss",&jwj_PTMiss_,"jwj_PTMiss/F");
    tree.Branch("jwj_PTMiss_trimmed1",&jwj_PTMiss_trimmed1_,"jwj_PTMiss_trimmed1/F");

}

float computeHT( std::vector< fastjet::PseudoJet > psjets, float cutoff = 0. ){
    
    float HT = 0;
    for (unsigned int i = 0; i < psjets.size(); i++){
        if (psjets[i].pt() > cutoff){
            HT += psjets[i].pt();
        }
    }
    return HT;
    
}

float computePTMiss( std::vector< fastjet::PseudoJet > psjets, float cutoff = 0. ){
    
    fastjet::PseudoJet ptvec(0,0,0,0);
    for (unsigned int i = 0; i < psjets.size(); i++){
        if (psjets[i].pt() > cutoff){
            ptvec = ptvec + psjets[i];
        }
    }
    return ptvec.pt();
    
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void eventAnalysis() {
    
    setTDRStyle();
    
    TFile fout_("output/outtree_20.root", "RECREATE");
    TTree* tree_gen_ = new TTree("tree_gen_", "tree_gen_");
    TTree* tree_pfchs_ = new TTree("tree_pfchs_", "tree_pfchs_");
    TTree* tree_pf_ = new TTree("tree_pf_", "tree_pf_");
    
    addBranches(*tree_gen_);
    addBranches(*tree_pfchs_);
    addBranches(*tree_pf_);
        
    TTree* tree_particles_ = new TTree("tree_particles_", "tree_particles_");    
    tree_particles_->Branch("pt",&pt_part,"pt/F");
    tree_particles_->Branch("eta",&eta_part,"eta/F");    

    initVars();
    
    int nEvts = 0;
    int maxEvents = 11;
    
    // fill vector of pseudojets
    std::vector < fastjet::PseudoJet > genParticles;
    std::vector < fastjet::PseudoJet > pfCHS;        
    std::vector < fastjet::PseudoJet > pf;    

    std::string filenameT = "/uscms_data/d2/ntran/physics/Jets/PUPPI/fromDaniele/samples/Zj_20.dat";
    std::cout << "Processing " << filenameT << std::endl;
    fin.open(filenameT.c_str());

    while(true){

        readEvent( genParticles, pfCHS, pf );

        if (nEvts > 0){
            
            initVars();
            analyzeEvent( genParticles, *tree_gen_ );
            initVars();
            analyzeEvent( pfCHS, *tree_pfchs_ );                
            initVars();
            analyzeEvent( pf, *tree_pf_ );   
            
            // plot event
            if (nEvts == 1 || nEvts == 10){
                char canvname[150];
                sprintf( canvname, "displays/dis_gen_%i", nEvts );
                plotEvent(genParticles, canvname);
                sprintf( canvname, "displays/dis_pfchs_%i", nEvts );
                plotEvent(pfCHS, canvname);
                sprintf( canvname, "displays/dis_pf_%i", nEvts );                
                plotEvent(pf, canvname);                
            }
                        
            genParticles.clear();
            pfCHS.clear();
            pf.clear();                

        }

        nEvts++;
        if (nEvts % 5 == 0) std::cout << "event no. = " << nEvts << std::endl;
        if (nEvts == maxEvents){ break; }
        
        if(fin.eof()) break;

    }

    fout_.cd();
    tree_gen_->Write();
    tree_pfchs_->Write();
    tree_pf_->Write();
    tree_particles_->Write();
    fout_.Close();
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void readEvent( std::vector< fastjet::PseudoJet > &gen, std::vector< fastjet::PseudoJet > &pfchs, std::vector< fastjet::PseudoJet > &pf ){
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;
    
    while(true){
        
        fin  >> npart >> px >> py >> pz >> e >> pdgid >> isCh >> isPU;         
        
        if (px == 0 && py == 0 && pz == 0 && e == 0){
            return;
        }        
        
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( px, py, pz, e );
        if (fabs(curPseudoJet.eta()) < 5){
            if (isPU == 0){
                gen.push_back( curPseudoJet );            
            }
            if ((isPU == 0) || (isPU == 1 && isCh == 0 && fabs(curPseudoJet.eta()) < 2.5) || (isPU == 1 && fabs(curPseudoJet.eta()) > 2.5)){
                pfchs.push_back( curPseudoJet );
            }
            pf.push_back( curPseudoJet );        
        }
        
        if(fin.eof()) break;
    }
    
}

/////////////////////////////////////////////////////////////////////

void analyzeEvent( std::vector < fastjet::PseudoJet > constits, TTree &tree ){
    
    // recluster on the fly....
    double rParam = 0.7;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);    
    
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    fastjet::AreaDefinition fjAreaDefinition_wGhosts( fastjet::active_area_explicit_ghosts, fjActiveArea );
    
    fastjet::ClusterSequenceArea* thisClustering_ = new fastjet::ClusterSequenceArea(constits, jetDef, fjAreaDefinition);
    fastjet::ClusterSequenceArea* thisClustering_wGhosts_ = new fastjet::ClusterSequenceArea(constits, jetDef, fjAreaDefinition_wGhosts);
//    fastjet::ClusterSequence* thisClustering_basic_ = new fastjet::ClusterSequence(constits, jetDef);
    
    std::vector<fastjet::PseudoJet> out_jets_ = sorted_by_pt(thisClustering_->inclusive_jets(5.0));
    std::vector<fastjet::PseudoJet> out_jets_wGhosts_ = sorted_by_pt(thisClustering_wGhosts_->inclusive_jets(5.0));    
//    out_jets_basic_ = sorted_by_pt(thisClustering_basic_->inclusive_jets(5.0));
    
    fastjet::Pruner pruner1( fastjet::cambridge_algorithm, 0.1, 0.5 );
    fastjet::Filter trimmer1( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03)) );
    
    // ------------------------
    // compute clustered and unclustered HT (sum of pT), pT_miss for trimmed and untrimmed jets
    std::vector<fastjet::PseudoJet> trimmedJets;
    for (unsigned int a = 0; a < out_jets_wGhosts_.size(); a++){
        trimmedJets.push_back( trimmer1( out_jets_wGhosts_[a] ) );
    }
    
    float jet_HT = computeHT( out_jets_, 5.0 );
    float jet_PTMiss = computePTMiss( out_jets_, 5.0 );
    float jet_HT_trimmed1 = computeHT( trimmedJets, 5.0 );
    float jet_PTMiss_trimmed1 = computePTMiss( trimmedJets, 5.0 );
    float prt_HT = computeHT( constits );
    float prt_PTMiss = computePTMiss( constits );
    
    // ------------------------
    // jwj quantities

    // Jet parameters. 
    double Rjet = rParam;
    double pTcut = 5;
    
    // Subjet parameters. Need for trimming.
    double Rsub = 0.2;
    double fcut = 0.03;
    
    Selector trimmer=SelectorEventTrimmer(Rjet,pTcut,Rsub,fcut);
    JetMultiplicity Nj(Rjet,pTcut), Nj_trim(Rjet,pTcut,Rsub,fcut);
    TransverseEnergy HT(Rjet,pTcut), HT_trim(Rjet,pTcut,Rsub,fcut);
    MissingTransverseEnergy HTmiss(Rjet,pTcut), HTmiss_trim(Rjet,pTcut,Rsub,fcut);
    
    float jwj_Nj = Nj(constits);
    float jwj_Nj_trimmed1 = Nj(constits);
    float jwj_HT = HT(constits);
    float jwj_HT_trimmed1 = HT_trim(constits);
    float jwj_PTMiss = HTmiss(constits);
    float jwj_PTMiss_trimmed1 = HTmiss_trim(constits);
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    // FILL IN THE TREE
    
    for (unsigned int i = 0; i < out_jets_.size(); i++){
        //std::cout << "jet " << i << ": " << out_jets_[i].m() << ", pt = " << out_jets_[i].pt() << ", eta = " << out_jets_[i].eta() << ", n constits = " << out_jets_[i].constituents().size() << std::endl;
        if (i < 10){
            jet_m_[i] = out_jets_[i].m();
            jet_pt_[i] = out_jets_[i].pt();            
            jet_area_[i] = out_jets_[i].area();            

            fastjet::PseudoJet prunedJet1 = pruner1( out_jets_wGhosts_[i] );
            fastjet::PseudoJet trimmedJet1 = trimmer1( out_jets_wGhosts_[i] );

            jet_m_pruned1_[i] = prunedJet1.m();
            jet_pt_pruned1_[i] = prunedJet1.pt();            
            jet_area_pruned1_[i] = prunedJet1.area();            
            
            jet_m_trimmed1_[i] = trimmedJet1.m();
            jet_pt_trimmed1_[i] = trimmedJet1.pt();            
            jet_area_trimmed1_[i] = trimmedJet1.area();            
            
            jet_nconstituents_[i] = out_jets_[i].constituents().size();            
        }
        else break;
    }
    
    // event quantities
    jet_njets_ = out_jets_.size();

    jet_HT_ = jet_HT;
    jet_PTMiss_ = jet_PTMiss;
    jet_HT_trimmed1_ = jet_HT_trimmed1;
    jet_PTMiss_trimmed1_ = jet_PTMiss_trimmed1;
    prt_HT_ = prt_HT;
    prt_PTMiss_ = prt_PTMiss;

    jwj_Nj_ = jwj_Nj;
    jwj_Nj_trimmed1_ = jwj_Nj_trimmed1;
    jwj_HT_ = jwj_HT;
    jwj_HT_trimmed1_ = jwj_HT_trimmed1;
    jwj_PTMiss_ = jwj_PTMiss;
    jwj_PTMiss_trimmed1_ = jwj_PTMiss_trimmed1;
    
    tree.Fill();
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    thisClustering_->delete_self_when_unused();
    thisClustering_wGhosts_->delete_self_when_unused();
//    thisClustering_basic_->delete_self_when_unused();
}

void plotEvent( std::vector < fastjet::PseudoJet > constits, char* name ){
    TH2F* h2d = new TH2F( "h2d",";#eta;#phi;e",100, -5,5, 100,0,2*TMath::Pi() );
    for (unsigned int i = 0; i < constits.size(); i++){
        int curBin = h2d->FindBin(constits[i].eta(),constits[i].phi());
        h2d->SetBinContent( curBin, h2d->GetBinContent(curBin) + constits[i].e() );
    }
    TCanvas* can = new TCanvas("can","can",800,800);
    h2d->Draw("BOX");
    char oname[96];
    sprintf(oname,"%s.eps",name);
    can->SaveAs(oname);
    delete h2d;
    delete can;
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
tdrStyle->SetPadLeftMargin(0.14);
tdrStyle->SetPadRightMargin(0.04);

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

tdrStyle->cd();

}

