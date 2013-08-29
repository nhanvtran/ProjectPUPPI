#! /usr/bin/env python
import os
import glob
import math
import array
import sys
import time

import ROOT

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadLeftMargin(0.16);


############################################################

if __name__ == '__main__':

    file = ROOT.TFile("output/outtree_20.root");
    tree_gen = file.Get("tree_gen_");
    tree_pfchs = file.Get("tree_pfchs_");
    tree_pf = file.Get("tree_pf_");    

    h_pt_gen = ROOT.TH1F("h_pt_gen",";p_{T} (GeV); count",25,0,1000);
    h_pt_pfchs = ROOT.TH1F("h_pt_pfchs",";p_{T} (GeV); count",25,0,1000);
    h_pt_pf = ROOT.TH1F("h_pt_pf",";p_{T} (GeV); count",25,0,1000);        

    h_m_gen = ROOT.TH1F("h_m_gen",";m (GeV); count",30,0,150);
    h_m_pfchs = ROOT.TH1F("h_m_pfchs",";m (GeV); count",30,0,150);
    h_m_pf = ROOT.TH1F("h_m_pf",";m (GeV); count",30,0,150);   

    h_nconst_gen = ROOT.TH1F("h_nconst_gen",";n constituents; count",30,0,150);
    h_nconst_pfchs = ROOT.TH1F("h_nconst_pfchs",";n constituents; count",30,0,150);
    h_nconst_pf = ROOT.TH1F("h_nconst_pf",";n constituents; count",30,0,150);   
    
    for i in range(tree_gen.GetEntries()):
        
        tree_gen.GetEntry(i);
        tree_pfchs.GetEntry(i);
        tree_pf.GetEntry(i);                
    
        h_pt_gen.Fill( tree_gen.pt[0] );
        h_pt_pfchs.Fill( tree_pfchs.pt[0] );
        h_pt_pf.Fill( tree_pf.pt[0] );                    

        h_m_gen.Fill( tree_gen.m[0] );
        h_m_pfchs.Fill( tree_pfchs.m[0] );
        h_m_pf.Fill( tree_pf.m[0] );                    

        h_nconst_gen.Fill( tree_gen.nconstituents[0] );
        h_nconst_pfchs.Fill( tree_pfchs.nconstituents[0] );
        h_nconst_pf.Fill( tree_pf.nconstituents[0] );                    
    
    h_pt_pfchs.SetLineColor( 2 );
    h_m_pfchs.SetLineColor( 2 );    
    h_nconst_pfchs.SetLineColor( 2 );    

    h_pt_pf.SetLineColor( 4 );
    h_m_pf.SetLineColor( 4 );    
    h_nconst_pf.SetLineColor( 4 );    
        
    # now plot

    can_m = ROOT.TCanvas("can_m","can_m",800,800);
    h_m_gen.Draw();
    h_m_pfchs.Draw("sames");
    h_m_pf.Draw("sames");        
    can_m.SaveAs("plots/can_m.eps");
        
    can_pt = ROOT.TCanvas("can_pt","can_pt",800,800);
    h_pt_gen.Draw();
    h_pt_pfchs.Draw("sames");
    h_pt_pf.Draw("sames");        
    can_pt.SaveAs("plots/can_pt.eps");            

    can_nconst = ROOT.TCanvas("can_nconst","can_nconst",800,800);
    h_nconst_gen.Draw();
    h_nconst_pfchs.Draw("sames");
    h_nconst_pf.Draw("sames");        
    can_nconst.SaveAs("plots/can_nconst.eps");            

    ##----------------------------------------------------
    ##----------------------------------------------------

#    tree_particles = file.Get("tree_particles_");    
#    h_particles_eta = ROOT.TH1F("h_particles_eta","; eta; count",30, -6,6 );
#    h_particles_pt = ROOT.TH1F("h_particles_pt","; p_{T} (GeV); count",20, 0,20 );
#
#    print "tree_particles.GetEntries() = ",tree_particles.GetEntries()
#    for i in range(tree_particles.GetEntries()):
#        if i%100000 == 0: print "i",i
#        tree_particles.GetEntry(i);
#        h_particles_eta.Fill( tree_particles.eta );
#        h_particles_pt.Fill( tree_particles.pt );        
#
#    cpart_eta = ROOT.TCanvas("cpart_eta","cpart_eta",800,800);
#    h_particles_eta.Draw();
#    cpart_eta.SaveAs("plots/cpart_eta.eps");
#        
#    cpart_pt = ROOT.TCanvas("cpart_pt","cpart_pt",800,800);
#    h_particles_pt.Draw();
#    cpart_pt.SaveAs("plots/cpart_pt.eps");



        