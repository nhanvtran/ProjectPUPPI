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

    file = ROOT.TFile("../output/outtree_80.root");
    tree_gen = file.Get("tree_gen_");
    tree_pfchs = file.Get("tree_pfchs_");
    tree_pf = file.Get("tree_pf_");    

    
    
    vars = ["m","m_trimmed1","m_pruned1","pt","pt_trimmed1","pt_pruned1","area","area_trimmed1","area_pruned1","nconstituents"]
    vars_x = ["m","mtrimmed1","mpruned1","pt","pt_trimmed1","pt_pruned1","area","area_trimmed1","area_pruned1","nconstituents"]
    vars_nbin = [20]*10;
    vars_lo = [   0,   0,   0,   0,   0,   0,    0,   0,  0,    0 ];
    vars_hi = [ 200, 200, 200, 800, 800, 800,    2,   2,  2,  200 ];
    
    h_gen = [];
    h_pfchs = [];
    h_pf = [];
    
    for i in range(len(vars)):
        h_gen.append( ROOT.TH1F("hgen_"+vars[i],";"+vars_x[i]+";n.d.",vars_nbin[i],vars_lo[i],vars_hi[i]) );
        h_pfchs.append( ROOT.TH1F("hpfchs_"+vars[i],";"+vars_x[i]+";n.d.",vars_nbin[i],vars_lo[i],vars_hi[i]) );   
        h_pf.append( ROOT.TH1F("hpf_"+vars[i],";"+vars_x[i]+";n.d.",vars_nbin[i],vars_lo[i],vars_hi[i]) );   
    
    for i in range(tree_gen.GetEntries()):
        
        tree_gen.GetEntry(i);
        tree_pfchs.GetEntry(i);
        tree_pf.GetEntry(i);        
        for j in range(len(vars)):
            h_gen[j].Fill( getattr( tree_gen, vars[j] )[0] );             
            h_pfchs[j].Fill( getattr( tree_pfchs, vars[j] )[0] );             
            h_pf[j].Fill( getattr( tree_pf, vars[j] )[0] );             
   
        
    # now plot

    ##### plot
    # scale and color
    for j in range(len(vars)):
        if h_gen[j].Integral() > 0:   h_gen[j].Scale( 1./h_gen[j].Integral() );
        if h_pfchs[j].Integral() > 0: h_pfchs[j].Scale( 1./h_pfchs[j].Integral() );
        if h_pf[j].Integral() > 0:    h_pf[j].Scale( 1./h_pf[j].Integral() );
        
        h_pfchs[j].SetLineColor( 2 );
        h_pf[j].SetLineColor( 4 );

    # on canvas
    for j in range(len(vars)):
        
        tmpcan = ROOT.TCanvas("can"+str(j),"can"+str(j),800,800);
        h_gen[j].SetMaximum( 1.2*max( h_gen[j].GetMaximum(), h_pfchs[j].GetMaximum(), h_pf[j].GetMaximum() ) );
        h_gen[j].SetMinimum( 0 );        
        h_gen[j].Draw();
        h_pfchs[j].Draw("sames");
        h_pf[j].Draw("sames");        
        tmpcan.SaveAs("figs/"+vars[j]+".eps");
        tmpcan.SaveAs("figs/"+vars[j]+".png");        

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



        