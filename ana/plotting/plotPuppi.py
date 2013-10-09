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
    tree_gen = file.Get("tree_gen");
    tree_pf = file.Get("tree_pf");
    tree_pf_tr = file.Get("tree_pf_tr");    
    tree_pf_cl = file.Get("tree_pf_cl");    
    
    types = ["gen","pf","pf_tr","pf_cl"];
    
    h_nJets = [];
    h_jet_mass = [];    
    h_jet_pt = [];        
    for i in range(4):
        h_nJets.append(ROOT.TH1F("h_nJets_"+types[i],"; n jets ;n.d.",10,0,10));
        h_jet_mass.append(ROOT.TH1F("h_jet_mass_"+types[i],"; mass ;n.d.",20,0,200));
        h_jet_pt.append(ROOT.TH1F("h_jet_pt"+types[i],"; pT ;n.d.",20,0,800));         

    trees = [tree_gen, tree_pf, tree_pf_tr, tree_pf_cl];

    for j in range(4):
        
        for i in range(tree_gen.GetEntries()):
            trees[j].GetEntry(i);
            h_nJets[j].Fill( trees[j].njets );
            h_jet_mass[j].Fill( trees[j].v_jet_m[0] );
            h_jet_pt[j].Fill( trees[j].v_jet_pt[0] );                    
            
            
    for j in range (4):
        h_nJets[j].Scale(1./h_nJets[j].Integral());
        h_jet_mass[j].Scale(1./h_jet_mass[j].Integral());
        h_jet_pt[j].Scale(1./h_jet_pt[j].Integral());

        h_nJets[j].SetLineColor(j+1);
        h_jet_mass[j].SetLineColor(j+1);
        h_jet_pt[j].SetLineColor(j+1);
        
        
    # plot
    can_nJets = ROOT.TCanvas("can_nJets","can_nJets",800,800);
    h_nJets[0].Draw();
    for j in range (1,4): h_nJets[j].Draw("sames");
    can_nJets.SaveAs("figs/can_nJets.eps");

    can_jet_mass = ROOT.TCanvas("can_jet_mass","can_jet_mass",800,800);
    h_jet_mass[0].Draw();
    for j in range (1,4): h_jet_mass[j].Draw("sames");
    can_jet_mass.SaveAs("figs/can_jet_mass.eps");

    can_jet_pt = ROOT.TCanvas("can_jet_pt","can_jet_pt",800,800);
    h_jet_pt[0].Draw();
    for j in range (1,4): h_jet_pt[j].Draw("sames");
    can_jet_pt.SaveAs("figs/can_jet_pt.eps");

                
#    
#    # jet vars
#    vars = ["jet_m","jet_m_trimmed1","jet_m_pruned1",
#            "jet_pt","jet_pt_trimmed1","jet_pt_pruned1",
#            "jet_area","jet_area_trimmed1","jet_area_pruned1",
#            "jet_nconstituents"]
#    vars_x = ["m","mtrimmed1","mpruned1","pt","pt_trimmed1","pt_pruned1","area","area_trimmed1","area_pruned1","nconstituents"]
#    vars_nbin = [20]*10;
#    vars_lo = [   0,   0,   0,   0,   0,   0,    0,   0,  0,    0 ];
#    vars_hi = [ 200, 200, 200, 800, 800, 800,    2,   2,  2,  200 ];
#    # event vars
#    evars = ["jet_njets",
#             "jet_HT", "jet_PTMiss","jet_HT_trimmed1", "jet_PTMiss_trimmed1",
#             "prt_HT","prt_PTMiss",
#             "jwj_Nj","jwj_Nj_trimmed1",
#             "jwj_HT","jwj_HT_trimmed1","jwj_PTMiss","jwj_PTMiss_trimmed1"];
#    evars_x = ["jet_njets",
#              "jet_HT", "jet_PTMiss","jet_HT_trimmed1", "jet_PTMiss_trimmed1",
#              "prt_HT","prt_PTMiss",
#              "jwj_Nj","jwj_Nj_trimmed1",
#              "jwj_HT","jwj_HT_trimmed1","jwj_PTMiss","jwj_PTMiss_trimmed1"]
#    evars_nbin = [20]*13;
#    evars_lo = [0,0,0,0,0,
#                0,0,0,0,0,
#                0,0,0];
#    evars_hi = [40, 
#                1500, 1500, 1500, 1500, 
#                1500, 1500,
#                  40,   40,
#                1500, 1500, 1500, 1500];
#    
#    h_gen = [];
#    h_pfchs = [];
#    h_pf = [];
#    
#    for i in range(len(vars)):
#        h_gen.append( ROOT.TH1F("hgen_"+vars[i],";"+vars_x[i]+";n.d.",vars_nbin[i],vars_lo[i],vars_hi[i]) );
#        h_pfchs.append( ROOT.TH1F("hpfchs_"+vars[i],";"+vars_x[i]+";n.d.",vars_nbin[i],vars_lo[i],vars_hi[i]) );   
#        h_pf.append( ROOT.TH1F("hpf_"+vars[i],";"+vars_x[i]+";n.d.",vars_nbin[i],vars_lo[i],vars_hi[i]) );   
#    for i in range(len(evars)):
#        h_gen.append( ROOT.TH1F("hgen_"+evars[i],";"+evars_x[i]+";n.d.",evars_nbin[i],evars_lo[i],evars_hi[i]) );
#        h_pfchs.append( ROOT.TH1F("hpfchs_"+evars[i],";"+evars_x[i]+";n.d.",evars_nbin[i],evars_lo[i],evars_hi[i]) );   
#        h_pf.append( ROOT.TH1F("hpf_"+evars[i],";"+evars_x[i]+";n.d.",evars_nbin[i],evars_lo[i],evars_hi[i]) );   
#    
#    for i in range(tree_gen.GetEntries()):
#        
#        tree_gen.GetEntry(i);
#        tree_pfchs.GetEntry(i);
#        tree_pf.GetEntry(i);        
#        for j in range(len(vars)):
#            h_gen[j].Fill( getattr( tree_gen, vars[j] )[0] );             
#            h_pfchs[j].Fill( getattr( tree_pfchs, vars[j] )[0] );             
#            h_pf[j].Fill( getattr( tree_pf, vars[j] )[0] );             
#        for j in range(len(vars),len(vars)+len(evars)):
#            kk = j - len(vars);
#            h_gen[j].Fill( getattr( tree_gen, evars[kk] ) );             
#            h_pfchs[j].Fill( getattr( tree_pfchs, evars[kk] ) );             
#            h_pf[j].Fill( getattr( tree_pf, evars[kk] ) );             
#   
#        
#    # now plot
#
#    ##### plot
#    # scale and color
#    for j in range(len(vars)+len(evars)):
#        if h_gen[j].Integral() > 0:   h_gen[j].Scale( 1./h_gen[j].Integral() );
#        if h_pfchs[j].Integral() > 0: h_pfchs[j].Scale( 1./h_pfchs[j].Integral() );
#        if h_pf[j].Integral() > 0:    h_pf[j].Scale( 1./h_pf[j].Integral() );
#        
#        h_pfchs[j].SetLineColor( 2 );
#        h_pf[j].SetLineColor( 4 );
#
#    # on canvas
#    for j in range(len(vars)+len(evars)):
#        
#        tmpcan = ROOT.TCanvas("can"+str(j),"can"+str(j),800,800);
#        h_gen[j].SetMaximum( 1.2*max( h_gen[j].GetMaximum(), h_pfchs[j].GetMaximum(), h_pf[j].GetMaximum() ) );
#        h_gen[j].SetMinimum( 0 );        
#        h_gen[j].Draw();
#        h_pfchs[j].Draw("sames");
#        h_pf[j].Draw("sames");    
#        if j < len(vars):
#            tmpcan.SaveAs("figs/"+vars[j]+".eps");
#            tmpcan.SaveAs("figs/"+vars[j]+".png");        
#        else:
#            tmpcan.SaveAs("figs/"+evars[j-len(vars)]+".eps");
#            tmpcan.SaveAs("figs/"+evars[j-len(vars)]+".png");

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



        