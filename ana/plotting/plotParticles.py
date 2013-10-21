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

    file = ROOT.TFile("../output/outtree_20.root");
    tree_particles = file.Get("tree_particles");

    h_puppiW_PU = ROOT.TH1F("h_puppiW_PU",";puppi weight; N",20,0,1);
    h_puppiW_LV = ROOT.TH1F("h_puppiW_LV",";puppi weight; N",20,0,1);    
    h_cleansedW_PU = ROOT.TH1F("h_cleansedW_PU",";cleansed weight; N",20,0,1);
    h_cleansedW_LV = ROOT.TH1F("h_cleansedW_LV",";cleansed weight; N",20,0,1);    

    tree_particles.Print();

    for i in range(tree_particles.GetEntries()):
        
        tree_particles.GetEntry(i);
        if tree_particles.p_isPU == 1 and tree_particles.p_puppiW >= 0:
            h_puppiW_PU.Fill( tree_particles.p_puppiW );
        if tree_particles.p_isPU == 0 and tree_particles.p_puppiW >= 0:
            h_puppiW_LV.Fill( tree_particles.p_puppiW );
        if tree_particles.p_isPU == 1 and tree_particles.p_cleansedW >= 0:
            h_cleansedW_PU.Fill( tree_particles.p_cleansedW );
        if tree_particles.p_isPU == 0 and tree_particles.p_cleansedW >= 0:
            h_cleansedW_LV.Fill( tree_particles.p_cleansedW );
        
    
    h_puppiW_PU.SetLineColor( 2 );
    h_puppiW_LV.SetLineColor( 1 );    
    h_puppiW_PU.SetLineWidth( 2 );
    h_puppiW_LV.SetLineWidth( 2 );    
    h_puppiW_PU.Scale(1./h_puppiW_PU.Integral());
    h_puppiW_LV.Scale(1./h_puppiW_LV.Integral());

    h_cleansedW_PU.SetLineColor( 2 );
    h_cleansedW_LV.SetLineColor( 1 );    
    h_cleansedW_PU.SetLineWidth( 2 );
    h_cleansedW_LV.SetLineWidth( 2 );    
    h_cleansedW_PU.SetLineStyle( 2 );
    h_cleansedW_LV.SetLineStyle( 2 );    
    h_cleansedW_PU.Scale(1./h_cleansedW_PU.Integral());
    h_cleansedW_LV.Scale(1./h_cleansedW_LV.Integral());
        
    canW = ROOT.TCanvas("canW","canW",800,800);
    h_puppiW_LV.SetMaximum( 1.2*max( h_puppiW_LV.GetMaximum(),h_puppiW_PU.GetMaximum(),h_cleansedW_LV.GetMaximum(),h_cleansedW_PU.GetMaximum() ) );
    h_puppiW_LV.Draw();
    h_puppiW_PU.Draw("sames");
    h_cleansedW_LV.Draw("sames");
    h_cleansedW_PU.Draw("sames");
    ROOT.gPad.SetLogy();    
    canW.SaveAs("figs_part/p_puppiW.eps");
    
    
    
    
    