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

def getmax( hs ):
    max = -999;
    for i in range(len(hs)):
        if hs[i] .GetMaximum() > max: max = hs[i].GetMaximum();
    return max;

if __name__ == '__main__':


    file = ROOT.TFile("../output/outtree_80.root");
    tree_gen = file.Get("tree_gen");
    tree_pf = file.Get("tree_pf");
    tree_pfchs = file.Get("tree_pfchs");
    tree_pf_tr = file.Get("tree_pf_tr");    
    tree_pf_cl = file.Get("tree_pf_cl");    
    tree_pf_puppi = file.Get("tree_pf_puppi");    
    
    types = ["gen","pf","pfchs","pf_tr","pf_cl","pf_puppi"];
    trees = [tree_gen, tree_pf, tree_pfchs, tree_pf_tr, tree_pf_cl,tree_pf_puppi];
    colors = [1,2,4,6,7,2];
    linestyles = [1,1,1,1,1,2];    
    nTypes = len(types);
    
    print "tree_gen.GetEntries() = ",tree_gen.GetEntries();
    print "tree_pf.GetEntries() = ",tree_pf.GetEntries();
    print "tree_pfchs.GetEntries() = ",tree_pfchs.GetEntries();
    print "tree_pf_tr.GetEntries() = ",tree_pf_tr.GetEntries();
    print "tree_pf_cl.GetEntries() = ",tree_pf_cl.GetEntries();
    print "tree_pf_puppi.GetEntries() = ",tree_pf_puppi.GetEntries();

    ## -----------------------------
    ## declare histograms
    h_nJets = [];
    h_jet_mass = [];    
    h_jet_pt = [];
    h_mass_corr = [];    
    for i in range(nTypes):
        h_nJets.append(ROOT.TH1F("h_nJets_"+types[i],"; n jets ;n.d.",50,0,50));
        h_jet_mass.append(ROOT.TH1F("h_jet_mass_"+types[i],"; mass ;n.d.",20,0,200));
        h_jet_pt.append(ROOT.TH1F("h_jet_pt_"+types[i],"; pT ;n.d.",20,0,800));         
        h_mass_corr.append(ROOT.TH2F("h_mass_corr_"+types[i],"; true mass; reco mass",100,0,200,100,0,200));

    ## -----------------------------
    ## fill histograms    
    for j in range(nTypes):
        
        for i in range(tree_gen.GetEntries()):
            trees[j].GetEntry(i);
            h_nJets[j].Fill( trees[j].njets );
            h_jet_mass[j].Fill( trees[j].v_jet_m[0] );
            h_jet_pt[j].Fill( trees[j].v_jet_pt[0] );        
            
    for j in range(trees[0].GetEntries()):
        for i in range(nTypes): trees[i].GetEntry(j);
        for i in range(nTypes):
            h_mass_corr[i].Fill(trees[0].v_jet_m[0],trees[i].v_jet_m[0]);
            
    ## -----------------------------
    ## properties of histograms            
    for j in range (nTypes):
        h_nJets[j].Scale(1./h_nJets[j].Integral());
        h_jet_mass[j].Scale(1./h_jet_mass[j].Integral());
        h_jet_pt[j].Scale(1./h_jet_pt[j].Integral());

        h_nJets[j].SetLineColor(colors[j]);
        h_jet_mass[j].SetLineColor(colors[j]);
        h_jet_pt[j].SetLineColor(colors[j]);
        h_nJets[j].SetLineWidth(2);
        h_jet_mass[j].SetLineWidth(2);
        h_jet_pt[j].SetLineWidth(2);
        
        h_nJets[j].SetLineStyle(linestyles[j]);
        h_jet_mass[j].SetLineStyle(linestyles[j]);
        h_jet_pt[j].SetLineStyle(linestyles[j]);
        
    leg = ROOT.TLegend(0.7,0.7,0.9,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    for j in range(nTypes):
        leg.AddEntry(h_nJets[j],types[j],"l");

    # plot
    can_nJets = ROOT.TCanvas("can_nJets","can_nJets",800,800);
    h_nJets[0].SetMaximum( 1.2*getmax(h_nJets) );
    h_nJets[0].Draw();
    for j in range (1,nTypes): h_nJets[j].Draw("sames");
    leg.Draw();
    can_nJets.SaveAs("figs/can_nJets.eps");

    can_jet_mass = ROOT.TCanvas("can_jet_mass","can_jet_mass",800,800);
    h_jet_mass[0].SetMaximum( 1.2*getmax(h_jet_mass) );
    h_jet_mass[0].Draw();
    for j in range (1,nTypes): h_jet_mass[j].Draw("sames");
    leg.Draw();    
    can_jet_mass.SaveAs("figs/can_jet_mass.eps");

    can_jet_pt = ROOT.TCanvas("can_jet_pt","can_jet_pt",800,800);
    h_jet_pt[0].SetMaximum( 1.2*getmax(h_jet_pt) );
    h_jet_pt[0].Draw();
    for j in range (1,nTypes): h_jet_pt[j].Draw("sames");
    leg.Draw();    
    can_jet_pt.SaveAs("figs/can_jet_pt.eps");

    for i in range(nTypes):
        print "corr factor ",i," = ", h_mass_corr[i].GetCorrelationFactor();   
        cantmp = ROOT.TCanvas("cancorr"+str(i),"cancorr"+str(i),800,800);
        h_mass_corr[i].Draw("box");
        ban1 = ROOT.TLatex(0.25,0.80,("r = "+str(h_mass_corr[i].GetCorrelationFactor()) ));
        ban1.SetNDC()
        ban1.SetTextSize(0.035)
        ban1.Draw();
        cantmp.SaveAs("figs/can_mcorr_"+types[i]+".eps");
        
        
        