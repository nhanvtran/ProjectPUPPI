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

// Global variables

// read in the file 
void setTDRStyle();

/////////////////////////////////////////////////////////////////////
// Tree variables
/////////////////////////////////////////////////////////////////////
ifstream fin;
int njets_;
std::vector<float> v_jet_m_;
std::vector<float> v_jet_pt_;
std::vector<float> v_jet_eta_;

/////////////////////////////////////////////////////////////////////
// Helper functions
/////////////////////////////////////////////////////////////////////
void readEvent( std::vector< fastjet::PseudoJet > &allParticles, std::vector<int> &v_isPU, std::vector<int> &v_isCh );
void analyzeEvent( std::vector < fastjet::PseudoJet > constits, TTree &tree, char* tag );
void plotEvent( std::vector < fastjet::PseudoJet > constits, char* name );


void initVars(){
    njets_ = -1.;
    v_jet_m_.clear();
    v_jet_pt_.clear();    
    v_jet_eta_.clear();        
}

void addBranches( TTree &tree ){
    tree.Branch("njets",&njets_);
    tree.Branch("v_jet_m",&v_jet_m_);
    tree.Branch("v_jet_pt",&v_jet_pt_);
    tree.Branch("v_jet_eta",&v_jet_eta_);    
}


/////////////////////////////////////////////////////////////////////
void puppiAnalyzer() {
        
    
    gROOT->ProcessLine("#include <vector>");    
    setTDRStyle();
            
    int nEvts = 0;
    int maxEvents = 1000;
    
    
    TFile fout_("output/outtree_20.root", "RECREATE");
    TTree* tree_gen = new TTree("tree_gen", "tree_gen");
    TTree* tree_pf = new TTree("tree_pf", "tree_pf");
    TTree* tree_pf_tr = new TTree("tree_pf_tr", "tree_pf_tr");
    TTree* tree_pf_cl = new TTree("tree_pf_cl", "tree_pf_cl");
    
    addBranches(*tree_gen);
    addBranches(*tree_pf);
    addBranches(*tree_pf_tr);    
    addBranches(*tree_pf_cl);        
    
    std::vector < fastjet::PseudoJet > allParticles;
    std::vector < int > v_isPU;
    std::vector < int > v_isCh;

    std::string filenameT = "/uscms_data/d2/ntran/physics/Jets/PUPPI/fromDaniele/samples/Zj_20.dat";
    std::cout << "Processing " << filenameT << std::endl;
    fin.open(filenameT.c_str());

    while(true){

        readEvent( allParticles, v_isPU, v_isCh );
        
        puppiContainer curEvent(allParticles, v_isPU, v_isCh);
        std::vector<fastjet::PseudoJet> trimmedParticles = curEvent.trimEvent();
        std::vector<fastjet::PseudoJet> cleansedParticles = curEvent.cleanseEvent();
        
        if (nEvts > 0){
            
            char canvname[150];
            
            if (nEvts < 10){
                sprintf( canvname, "displays/v2/dis_gen_%i", nEvts );
                plotEvent( curEvent.genParticles(), canvname );
                sprintf( canvname, "displays/v2/dis_pf_%i", nEvts );
                plotEvent( curEvent.pfParticles(), canvname );
                sprintf( canvname, "displays/v2/dis_pf_tr_%i", nEvts );
                plotEvent( trimmedParticles, canvname );
                sprintf( canvname, "displays/v2/dis_pf_cl_%i", nEvts );
                plotEvent( cleansedParticles, canvname );
            }
            
            sprintf( canvname, "dummy" );
            initVars();
            analyzeEvent( curEvent.genParticles(), *tree_gen, canvname );
            initVars();
            analyzeEvent( curEvent.pfParticles(), *tree_pf, canvname );
            initVars();
            analyzeEvent( trimmedParticles, *tree_pf_tr, canvname );
            initVars();
            analyzeEvent( cleansedParticles, *tree_pf_cl, canvname );
            
            allParticles.clear();
            v_isPU.clear();
            v_isCh.clear();
        }
            
        nEvts++;
        if (nEvts % 100 == 0) std::cout << "event no. = " << nEvts << std::endl;
        if (nEvts == maxEvents){ break; }

        if(fin.eof()) break;

    }

    fout_.cd();
    tree_gen->Write();
    tree_pf->Write();
    tree_pf_tr->Write();
    tree_pf_cl->Write();
    fout_.Close();
    
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
void readEvent( std::vector< fastjet::PseudoJet > &allParticles, std::vector<int> &v_isPU, std::vector<int> &v_isCh ){
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;
    
    while(true){
        
        fin  >> npart >> px >> py >> pz >> e >> pdgid >> isCh >> isPU;         
        
        if (px == 0 && py == 0 && pz == 0 && e == 0){
            return;
        }        
        
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( px, py, pz, e );
        if (fabs(curPseudoJet.eta()) < 5){
            allParticles.push_back( curPseudoJet );
            v_isPU.push_back(isPU);
            v_isCh.push_back(isCh);            
        }
        
        if(fin.eof()) break;
    }
    
}

void plotEvent( std::vector < fastjet::PseudoJet > constits, char* name ){
    
    TH2F* h2d = new TH2F( "h2d",";#eta;#phi;e",100, -5,5, 100,0,2*TMath::Pi() );
    for (unsigned int i = 0; i < constits.size(); i++){
        int curBin = h2d->FindBin(constits[i].eta(),constits[i].phi());
        h2d->SetBinContent( curBin, h2d->GetBinContent(curBin) + constits[i].e() );
    }
        
    TCanvas* can = new TCanvas("can","can",1000,800);
    h2d->Draw("BOX");

    char oname[96];
    sprintf(oname,"%s.eps",name);
    sprintf(oname,"%s.png",name);
    can->SaveAs(oname);
    delete h2d;
    delete can;
    
}

void analyzeEvent( std::vector < fastjet::PseudoJet > constits, TTree &tree, char* tag ){
    
    //std::cout << "constits.size() = " << constits.size() << std::endl;
    
    // recluster on the fly....
    double rParam = 0.7;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);    
    
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    fastjet::ClusterSequenceArea* thisClustering_ = new fastjet::ClusterSequenceArea(constits, jetDef, fjAreaDefinition);
    std::vector<fastjet::PseudoJet> out_jets_ = sorted_by_pt(thisClustering_->inclusive_jets(20.0));
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    // FILL IN THE TREE
    //std::cout << "out_jets_.size() = " << out_jets_.size() << std::endl; 
    njets_ = (int) out_jets_.size();
    for (unsigned int i = 0; i < out_jets_.size(); i++){
        v_jet_m_.push_back( out_jets_[i].m() );
        v_jet_pt_.push_back( out_jets_[i].pt() );        
        v_jet_eta_.push_back( out_jets_[i].eta() );                
    }
    
    // event quantities    
    tree.Fill();
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    thisClustering_->delete_self_when_unused();

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

