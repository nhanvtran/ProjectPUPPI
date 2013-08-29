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

#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

using namespace fastjet;


// mass quantities
float m_[10];
float pt_[10];

// Groomed quantities
//     trimmed1 = 
//     pruned1  = 
float m_trimmed1_[10];
float m_pruned1_[10];
float pt_trimmed1_[10];
float pt_pruned1_[10];
float area_[10];
float area_trimmed1_[10];
float area_pruned1_[10];
float nconstituents_[10];

//other
float njets_;

///////////////////per particle
float pt_part;
float eta_part;

void clusterEvent( std::vector < fastjet::PseudoJet > constits, std::vector< float > pdgids, TTree &tree );

void initVars(){
    for (int i = 0; i < 10; i++ ){
        m_[i] = -1.;
        pt_[i] = -1.;   
        nconstituents_[i] = -1.;
    }
    njets_ = -1.;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void eventAnalysis() {
    
    TFile fout_("output/outtree_80.root", "RECREATE");
    TTree* tree_gen_ = new TTree("tree_gen_", "tree_gen_");
    TTree* tree_pfchs_ = new TTree("tree_pfchs_", "tree_pfchs_");
    TTree* tree_pf_ = new TTree("tree_pf_", "tree_pf_");
    TTree* tree_particles_ = new TTree("tree_particles_", "tree_particles_");
    
    tree_gen_->Branch("m",&m_,"m[10]/F");
    tree_pfchs_->Branch("m",&m_,"m[10]/F");
    tree_pf_->Branch("m",&m_,"m[10]/F");
    tree_gen_->Branch("pt",&pt_,"pt[10]/F");
    tree_pfchs_->Branch("pt",&pt_,"pt[10]/F");
    tree_pf_->Branch("pt",&pt_,"pt[10]/F");
    tree_gen_->Branch("nconstituents",&nconstituents_,"nconstituents[10]/F");
    tree_pfchs_->Branch("nconstituents",&nconstituents_,"nconstituents[10]/F");
    tree_pf_->Branch("nconstituents",&nconstituents_,"nconstituents[10]/F");
    
    tree_gen_->Branch("njets",&njets_,"njets/F");
    tree_pfchs_->Branch("njets",&njets_,"njets/F");
    tree_pf_->Branch("njets",&njets_,"njets/F");
    
    tree_particles_->Branch("pt",&pt_part,"pt/F");
    tree_particles_->Branch("eta",&eta_part,"eta/F");    
    
    initVars();
    
    // read in the file 
    ifstream fin;
    std::string filenameT = "/uscms_data/d2/ntran/physics/Jets/PUPPI/fromDaniele/samples/Zj_80.dat";
    std::cout << "Processing " << filenameT << std::endl;
    fin.open(filenameT.c_str());

    int nEvts = 0;
    int maxEvents = 500;
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;
    
    // fill vector of pseudojets
    std::vector < fastjet::PseudoJet > genParticles;
    std::vector < fastjet::PseudoJet > pfCHS;        
    std::vector < fastjet::PseudoJet > pf;    
    std::vector < float > genParticles_pdgids;
    std::vector < float > pfCHS_pdgids;        
    std::vector < float > pf_pdgids;
    
    while (true) {
        
        fin  >> npart >> px >> py >> pz >> e >> pdgid >> isCh >> isPU;         
        
        if (px == 0 && py == 0 && pz == 0 && e == 0){

            // process the last event
            if (nEvts > 0){
//                std::cout << "------------------ event " << nEvts << std::endl;
//                std::cout << "genParticles.size() = " << genParticles.size() << std::endl;
                initVars();
                clusterEvent( genParticles, genParticles_pdgids, *tree_gen_ );
//                std::cout << "pfCHS.size() = " << pfCHS.size() << std::endl;
                initVars();
                clusterEvent( pfCHS, pfCHS_pdgids, *tree_pfchs_ );                
//                std::cout << "pf.size() = " << pf.size() << std::endl;    
                initVars();
                clusterEvent( pf, pf_pdgids, *tree_pf_ );      
                
                /* turn off
                for (unsigned int k = 0; k < pf.size(); k++){
                    pt_part = pf[k].pt();
                    eta_part = pf[k].eta();
                    tree_particles_->Fill();
                }
                */ 
                
                // delete vectors
                genParticles.clear();
                pfCHS.clear();
                pf.clear();
                genParticles_pdgids.clear();
                pfCHS_pdgids.clear();
                pf_pdgids.clear();                
            }
            
            nEvts++;
            if (nEvts % 25 == 0) std::cout << "event no. = " << nEvts << std::endl;
            if (nEvts == maxEvents){ break; }
            continue;
        }        
        
        if(fin.eof()) break;
        
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( px, py, pz, e );
        if (fabs(curPseudoJet.eta()) < 5){
            if (isPU == 0){
                genParticles.push_back( curPseudoJet );            
                genParticles_pdgids.push_back( pdgid );                        
            }
            if ((isPU == 0) || (isPU == 1 && isCh == 0 && fabs(curPseudoJet.eta()) < 2.5)){
                pfCHS.push_back( curPseudoJet );
                pfCHS_pdgids.push_back( pdgid );                        
            }
            pf.push_back( curPseudoJet );        
            pf_pdgids.push_back( pdgid );                                
        }
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

void clusterEvent( std::vector < fastjet::PseudoJet > constits, std::vector< float > pdgids, TTree &tree ){
    
    // recluster on the fly....
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, 0.7);    
    
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
    
//    std::cout << "size of outjets = " << out_jets_.size() << std::endl;

    fastjet::Pruner pruner1( fastjet::cambridge_algorithm, 0.1, 0.5 );
    fastjet::Filter trimmer1( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03)) );
    
    // FILL IN THE TREE
    
    for (unsigned int i = 0; i < out_jets_.size(); i++){
        //std::cout << "jet " << i << ": " << out_jets_[i].m() << ", pt = " << out_jets_[i].pt() << ", eta = " << out_jets_[i].eta() << ", n constits = " << out_jets_[i].constituents().size() << std::endl;
        if (i < 10){
            m_[i] = out_jets_[i].m();
            pt_[i] = out_jets_[i].pt();            
            area_[i] = out_jets_[i].area();            

            fastjet::PseudoJet prunedJet1 = pruner1( out_jets_wGhosts_[i] );
            fastjet::PseudoJet trimmedJet1 = trimmer1( out_jets_wGhosts_[i] );

            m_pruned1_[i] = prunedJet1.m();
            pt_pruned1_[i] = prunedJet1.pt();            
            area_pruned1_[i] = prunedJet1.area();            
            
            m_trimmed1_[i] = trimmedJet1.m();
            pt_trimmed1_[i] = trimmedJet1.pt();            
            area_trimmed1_[i] = trimmedJet1.area();            
            
            nconstituents_[i] = out_jets_[i].constituents().size();            
        }
        else break;
    }
    njets_ = out_jets_.size();

    tree.Fill();
    
    thisClustering_->delete_self_when_unused();
    thisClustering_wGhosts_->delete_self_when_unused();
//    thisClustering_basic_->delete_self_when_unused();
}

