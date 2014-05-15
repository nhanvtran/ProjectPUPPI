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
#include <fastjet/JetDefinition.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include <fastjet/tools/GridMedianBackgroundEstimator.hh>


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
int njets_corr_;

float rho_;
float sigma_;

float sumEt_;
float etMissX_;
float etMissY_;

std::vector<float> v_jet_m_;
std::vector<float> v_jet_pt_;
std::vector<float> v_jet_eta_;
std::vector<float> v_jet_phi_;

// trimmed info
std::vector<float> v_jet_m_trimmed_;
std::vector<float> v_jet_pt_trimmed_;
std::vector<float> v_jet_eta_trimmed_;
std::vector<float> v_jet_phi_trimmed_;

// 4-vector subtracted
std::vector<float> v_jet_m_4Vcorr_;
std::vector<float> v_jet_pt_4Vcorr_;
std::vector<float> v_jet_eta_4Vcorr_;
std::vector<float> v_jet_phi_4Vcorr_;

float p_isPU, p_isCH, p_px, p_py, p_pz, p_e, p_puppiW_pfchs, p_cleansedW, p_puppiW_chLV, p_puppiW_all;
float p_alphas_chLV, p_alphas_all;

// some randoms
double curRho;
double pfRho;
double pfchsRho;
bool isPFCHS;

/////////////////////////////////////////////////////////////////////
// Helper functions
/////////////////////////////////////////////////////////////////////
void readEvent( std::vector< fastjet::PseudoJet > &allParticles, std::vector<int> &v_isPU, std::vector<int> &v_isCh );
std::vector< fastjet::PseudoJet > analyzeEvent( std::vector < fastjet::PseudoJet > constits, TTree &tree, char* tag, double vRparam );
void plotEvent( std::vector < fastjet::PseudoJet > constits, char* name, std::vector < fastjet::PseudoJet > jets );


void initVars(){
    njets_ = -1.;
    njets_corr_ = -1.;
    
    rho_ = -1.;
    sigma_ = -1.;
    
    sumEt_ = -1.;
    etMissX_ = -1.;
    etMissY_ = -1.;    
    
    v_jet_m_.clear();
    v_jet_pt_.clear();    
    v_jet_eta_.clear(); 
    v_jet_phi_.clear(); 
    v_jet_m_trimmed_.clear();
    v_jet_pt_trimmed_.clear();
    v_jet_eta_trimmed_.clear();
    v_jet_phi_trimmed_.clear();
    v_jet_m_4Vcorr_.clear();
    v_jet_pt_4Vcorr_.clear();
    v_jet_eta_4Vcorr_.clear();
    v_jet_phi_4Vcorr_.clear();
}

void addBranches( TTree &tree ){

    tree.Branch("njets",&njets_);
    tree.Branch("njets_corr",&njets_corr_);

    tree.Branch("rho",&rho_);
    tree.Branch("sigma",&sigma_);
    
    tree.Branch("sumEt",&sumEt_);
    tree.Branch("etMissX",&etMissX_);
    tree.Branch("etMissY",&etMissY_);
    
    tree.Branch("v_jet_m",&v_jet_m_);
    tree.Branch("v_jet_pt",&v_jet_pt_);
    tree.Branch("v_jet_eta",&v_jet_eta_);  
    tree.Branch("v_jet_phi",&v_jet_phi_);  

    tree.Branch("v_jet_m_trimmed",&v_jet_m_trimmed_);    
    tree.Branch("v_jet_pt_trimmed",&v_jet_pt_trimmed_);
    tree.Branch("v_jet_eta_trimmed",&v_jet_eta_trimmed_);
    tree.Branch("v_jet_phi_trimmed",&v_jet_phi_trimmed_);
    
    tree.Branch("v_jet_m_4Vcorr",&v_jet_m_4Vcorr_);    
    tree.Branch("v_jet_pt_4Vcorr",&v_jet_pt_4Vcorr_);    
    tree.Branch("v_jet_eta_4Vcorr",&v_jet_eta_4Vcorr_);    
    tree.Branch("v_jet_phi_4Vcorr",&v_jet_phi_4Vcorr_);   
    
}

void initParticleVars(){
    p_isPU = -1.;
    p_isCH = -1.;
    
    p_px = -1.;
    p_py = -1.;
    p_pz = -1.;
    p_e = -1.;
    
    p_puppiW_pfchs = -1.;
    p_puppiW_chLV = -1.;
    p_puppiW_all = -1.;

    p_alphas_chLV = -1.;
    p_alphas_all = -1.;
    
    
    p_cleansedW = -1.;            
}

void addParticleBranches( TTree &tree ){
    tree.Branch("p_isPU",&p_isPU);
    tree.Branch("p_isCH",&p_isCH);
    tree.Branch("p_px",&p_px);
    tree.Branch("p_py",&p_py);
    tree.Branch("p_pz",&p_pz);
    tree.Branch("p_e",&p_e);
    tree.Branch("p_puppiW_pfchs",&p_puppiW_pfchs);
    tree.Branch("p_puppiW_chLV",&p_puppiW_chLV);
    tree.Branch("p_puppiW_all",&p_puppiW_all);
    tree.Branch("p_alphas_chLV",&p_alphas_chLV);
    tree.Branch("p_alphas_all",&p_alphas_all);
    tree.Branch("p_cleansedW",&p_cleansedW);
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/// m a i n   m a i n   m a i n   m a i n   m a i n   m a i n 
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

int main( int argc, char **argv ) {
    
    char filename[30];
    int maxEventsVal = atoi(argv[1]);
    strcpy( filename, argv[2] );
    int NPU = atoi(argv[3]);
    
    gROOT->ProcessLine("#include <vector>");    
    setTDRStyle();
    
    int nEvts = 0;
    int maxEvents = maxEventsVal;
    
    char tfname[150];
    sprintf( tfname, "output/outtree_%s_051514.root", filename );
    TFile fout_(tfname, "RECREATE");
    TTree* tree_gen = new TTree("tree_gen", "tree_gen");
    TTree* tree_pf = new TTree("tree_pf", "tree_pf");
    TTree* tree_pfchs = new TTree("tree_pfchs", "tree_pfchs");
    TTree* tree_pf_tr = new TTree("tree_pf_tr", "tree_pf_tr");
    TTree* tree_pf_cl = new TTree("tree_pf_cl", "tree_pf_cl");
    TTree* tree_pf_puppi = new TTree("tree_pf_puppi", "tree_pf_puppi");
    
    addBranches(*tree_gen);
    addBranches(*tree_pf);
    addBranches(*tree_pfchs);
    addBranches(*tree_pf_tr);    
    addBranches(*tree_pf_cl);   
    addBranches(*tree_pf_puppi);  
    
    TTree* tree_particles = new TTree("tree_particles", "tree_particles");
    addParticleBranches(*tree_particles);  
    
    std::vector < fastjet::PseudoJet > allParticles;
    std::vector < int > v_isPU;
    std::vector < int > v_isCh;
    
    char fname[150];
    sprintf( fname, "samples/%s.dat", filename );
    //sprintf( fname, "samples/minbias_%i.dat", PUscenario );
    std::cout << "Processing " << fname << std::endl;
    fin.open(fname);
    
    while(true){
        
        readEvent( allParticles, v_isPU, v_isCh );
        
        if (nEvts > 1){
            
            puppiContainer curEvent(allParticles, v_isPU, v_isCh, 1);
            
            std::vector<fastjet::PseudoJet> genParticles = curEvent.genParticles();
            std::vector<fastjet::PseudoJet> pfParticles = curEvent.pfParticles();
            std::vector<fastjet::PseudoJet> pfchsParticles = curEvent.pfchsParticles();
            //std::vector<fastjet::PseudoJet> trimmedParticles = curEvent.trimEvent();
            //std::vector<fastjet::PseudoJet> cleansedParticles = curEvent.cleanseEvent(0.3);
            //std::vector<fastjet::PseudoJet> puppiParticles = curEvent.puppiEvent(0.5,2.);
            
            std::vector<fastjet::PseudoJet> puppiParticles = curEvent.puppiEvent(NPU);
            
//            for (unsigned int kk = 0; kk < puppiParticles.size(); kk++){
//                std::cout << "puppiParticles.user_index() = " << puppiParticles[kk].user_index() << std::endl;
//            }
            
            char canvname[150];
                        
            //            sprintf( canvname, "dummy" );
            
            curRho   = -1;
            pfRho    = -1;
            pfchsRho = -1;
            //std::cout << "--------------------" << std::endl;
            initVars();
            //std::cout << "analyze gen" << std::endl;
            std::vector<fastjet::PseudoJet> genJets = analyzeEvent( genParticles, *tree_gen, canvname, 0.7 );
            initVars();
            //std::cout << "analyze pf" << std::endl;
            std::vector<fastjet::PseudoJet> pfJets = analyzeEvent( pfParticles, *tree_pf, canvname, 0.7 );
            pfRho = curRho;
            initVars();
            //std::cout << "analyze pfchs" << std::endl;
            isPFCHS = true;
            std::vector<fastjet::PseudoJet> pfchsJets = analyzeEvent( curEvent.pfchsParticles(), *tree_pfchs, canvname, 0.7 );
            pfchsRho = curRho;
            isPFCHS = false;
            //            initVars();
            //            analyzeEvent( trimmedParticles, *tree_pf_tr, canvname, 0.7 );
            //            initVars();
            //            analyzeEvent( cleansedParticles, *tree_pf_cl, canvname, 0.7 );
            initVars();
            //std::cout << "analyze puppi" << std::endl;
            std::vector<fastjet::PseudoJet> puppiJets = analyzeEvent( puppiParticles, *tree_pf_puppi, canvname, 0.7 );
            //
            
            
            if (nEvts < 200 and nEvts > 0){
                                                
                std::vector<float> puppiWeights_pfchs = curEvent.getPuppiWeights_pfchs();
                std::vector<float> puppiWeights_chLV = curEvent.getPuppiWeights_chLV();                
                std::vector<float> puppiWeights_all = curEvent.getPuppiWeights_all();                
                std::vector<float> alphas_chLV = curEvent.getPuppiAlphas_chLV();
                std::vector<float> alphas_all = curEvent.getPuppiAlphas_all();
                /////std::vector<float> cleansedWeights = curEvent.getCleansedWeights();
                // fill weights
                for (unsigned int a = 0; a < pfParticles.size(); a++){
                    if (pfParticles[a].user_index() == 1 || pfParticles[a].user_index() == 3) p_isPU = 1;
                    else p_isPU = 0;
                    if (pfParticles[a].user_index() == 2 || pfParticles[a].user_index() == 3) p_isCH = 1;
                    else p_isCH = 0;
                    p_px = pfParticles[a].px();
                    p_py = pfParticles[a].py();
                    p_pz = pfParticles[a].pz();
                    p_e = pfParticles[a].e();
//                    p_puppiW_pfchs = puppiWeights_pfchs[a];
                    p_puppiW_chLV = puppiWeights_chLV[a];
                    p_puppiW_all = puppiWeights_all[a];
                    p_alphas_chLV = alphas_chLV[a];
                    p_alphas_all = alphas_all[a];
                    /////p_cleansedW = cleansedWeights[a];
                    tree_particles->Fill();
                }
                
                if (nEvts < 0){
                                
                    sprintf( canvname, "displays/dis_gen_%i_%i", nEvts, NPU );
                    plotEvent( genParticles, canvname, genJets );
                    sprintf( canvname, "displays/dis_pf_%i_%i", nEvts, NPU );
                    plotEvent( pfParticles, canvname, pfJets);
                    sprintf( canvname, "displays/dis_pfchs_%i_%i", nEvts, NPU );
                    plotEvent( pfchsParticles, canvname, pfchsJets );                
                    //                sprintf( canvname, "displays/dis_pf_tr_%i_%i", nEvts, PUscenario );
                    //                plotEvent( trimmedParticles, canvname );
                    //                sprintf( canvname, "displays/dis_pf_cl_%i_%i", nEvts, PUscenario );
                    //                plotEvent( cleansedParticles, canvname );
                    sprintf( canvname, "displays/dis_pf_puppi_%i_%i", nEvts, NPU );
                    plotEvent( puppiParticles, canvname, puppiJets );
                }
            }            
            
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
    tree_pfchs->Write();    
    tree_pf_tr->Write();
    tree_pf_cl->Write();
    tree_pf_puppi->Write();
    tree_particles->Write();
    fout_.Close();
    
    return 0;
    
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
void readEvent( std::vector< fastjet::PseudoJet > &allParticles, std::vector<int> &v_isPU, std::vector<int> &v_isCh ){
    
    // hard-coded! 
    float etaMax = 5.;
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;
    
    while(true){
        
        fin  >> npart >> px >> py >> pz >> e >> pdgid >> isCh >> isPU;         
        
        if (px == 0 && py == 0 && pz == 0 && e == 0){
            return;
        }        
        
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( px, py, pz, e );
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
            v_isPU.push_back(isPU);
            v_isCh.push_back(isCh);            
        }
        
        if(fin.eof()) break;
    }
    
}

void plotEvent( std::vector < fastjet::PseudoJet > constits, char* name, std::vector < fastjet::PseudoJet > jets ){
    
    double maxCell = 0;
    
    TH2F* h2d = new TH2F( "h2d",";#eta;#phi;pT (GeV)",50, -5,5, 50,0,2*TMath::Pi() );
    TH2F* h2d_PU = new TH2F( "h2d_PU",";#eta;#phi;pT (GeV)",50, -5,5, 50,0,2*TMath::Pi() );
    
    for (unsigned int i = 0; i < constits.size(); i++){
        int curBin = h2d->FindBin(constits[i].eta(),constits[i].phi());
        if (constits[i].pt() > maxCell) maxCell = constits[i].pt();

        if (constits[i].user_index() == 1 || constits[i].user_index() == 3){
            h2d_PU->SetBinContent( curBin, h2d_PU->GetBinContent(curBin) + constits[i].pt() );
        }
        else{
            h2d->SetBinContent( curBin, h2d->GetBinContent(curBin) + constits[i].pt() );
        }
    }

    //std::cout << "maxCell = " << maxCell << std::endl;
    
    TCanvas* can = new TCanvas("can","can",1100,800);
    
    h2d->SetMaximum( 1.1*max(h2d->GetMaximum(),h2d_PU->GetMaximum()) );
    h2d->SetMinimum( 1e-1 );
//    h2d->SetLineWidth( 2 );
    h2d->Draw("COLZ");    
    h2d_PU->SetLineWidth( 1 );
    h2d_PU->SetLineColor( 14 );
    h2d_PU->Draw("BOX");
    h2d->Draw("COLZ SAME");
    can->SetLogz();

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
            cir->SetLineWidth( 2 );
        }
        
        cir->Draw("sames");
    }
    
    char oname[96];
    sprintf(oname,"%s.pdf",name);
    can->SaveAs(oname);
    sprintf(oname,"%s.png",name);
    can->SaveAs(oname);
    delete h2d;
    delete h2d_PU;
    delete can;
}

std::vector< fastjet::PseudoJet > analyzeEvent( std::vector < fastjet::PseudoJet > constits, TTree &tree, char* tag, double vRparam ){
    
    //std::cout << "constits.size() = " << constits.size() << std::endl;
    
    // recluster on the fly....
    double rParam = vRparam;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);    
    
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area_explicit_ghosts, fjActiveArea );
    fastjet::ClusterSequenceArea* thisClustering_ = new fastjet::ClusterSequenceArea(constits, jetDef, fjAreaDefinition);
    std::vector<fastjet::PseudoJet> out_jets_ = sorted_by_pt(thisClustering_->inclusive_jets(15.0));
    
    // trim jet
    fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.3), fastjet::SelectorPtFractionMin(0.05)));
    
    // 4-vector subtraction
    std::vector<fastjet::PseudoJet> neutrals;
    float etarange = 5.0;
    if(isPFCHS){
        // get the neutrals and compute rho only inside of 2.5
        etarange = 2.5;
        for (unsigned int aa = 0; aa < constits.size(); aa++){
            if(constits[aa].user_index() <= 1) neutrals.push_back(constits[aa]);
        }
    }
    fastjet::GridMedianBackgroundEstimator lGrid(etarange,0.8);
    if (isPFCHS) lGrid.set_particles(neutrals);
    else lGrid.set_particles(constits);
    
    curRho = lGrid.rho();
    rho_ = lGrid.rho();
    sigma_ = lGrid.sigma();

    double forwardRho = lGrid.rho();
    double centralRho = lGrid.rho();
    if(isPFCHS){ forwardRho = pfRho; }
    //std::cout << "centralRho = " << centralRho << ", forwardRho = " << forwardRho << std::endl;
    
    njets_corr_ = 0;
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    // FILL IN THE TREE
    //std::cout << "out_jets_.size() = " << out_jets_.size() << std::endl;
    std::vector<fastjet::PseudoJet> corrjets;
    
    njets_ = (int) out_jets_.size();
    for (unsigned int i = 0; i < out_jets_.size(); i++){

        v_jet_m_.push_back( out_jets_[i].m() );
        v_jet_pt_.push_back( out_jets_[i].pt() );        
        v_jet_eta_.push_back( out_jets_[i].eta() );                
        v_jet_phi_.push_back( out_jets_[i].phi() );                
        
        // 4-vector subtraction
        PseudoJet pCorrJet = out_jets_.at(i);        
        PseudoJet pArea = out_jets_.at(i).area_4vector();
    
        double tmpRho = -9999;
        if (fabs(pCorrJet.eta()) > 2.5) tmpRho = forwardRho;
        else tmpRho = centralRho;
        
        pCorrJet -= tmpRho * pArea;
    
        v_jet_m_4Vcorr_.push_back(pCorrJet.m());
        v_jet_pt_4Vcorr_.push_back(pCorrJet.pt());
        v_jet_eta_4Vcorr_.push_back(pCorrJet.eta());
        v_jet_phi_4Vcorr_.push_back(pCorrJet.phi());

        // trim jet + 4-vector subtraction
        fastjet::PseudoJet trimmedJet = (trimmer)(out_jets_.at(i));
        PseudoJet tjetArea = trimmedJet.area_4vector();
        trimmedJet -= tmpRho*tjetArea;
        v_jet_m_trimmed_.push_back( trimmedJet.m() );
        v_jet_pt_trimmed_.push_back( trimmedJet.pt() );
        v_jet_eta_trimmed_.push_back( trimmedJet.eta() );
        v_jet_phi_trimmed_.push_back( trimmedJet.phi() );
        
        if (v_jet_pt_4Vcorr_[i] > 25.0){
            njets_corr_++;
            corrjets.push_back(pCorrJet);
        }
    }
    
    // MET variables
    float sumEt = 0;
    PseudoJet METvec(0,0,0,0);
    for (unsigned int i = 0; i < constits.size(); i++){
        sumEt += constits[i].Et();
        METvec += constits[i];
    }
    
    sumEt_ = sumEt;
    etMissX_ = fabs(-METvec.px());
    etMissY_ = fabs(-METvec.py());
    
    // event quantities    
    tree.Fill();
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    if (out_jets_.size() > 0) thisClustering_->delete_self_when_unused();
    
    return corrjets;
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
