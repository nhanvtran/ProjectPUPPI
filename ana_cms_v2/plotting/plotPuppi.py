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

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--nPU',action="store",type="int",dest="nPU",default=22)

(options, args) = parser.parse_args()

############################################################
def makeCanvas(hists, names, canname, isLog=False, norm=False, wGen=False):
    
    directory = "plots";
    colors = [1,4,2,6,7];
    widths = [3,3,3,3,3];
    styles = [1,1,2,1,1];
    
    max = -999.;
    for hist in hists:
        if norm: hist.Scale(1./hist.Integral());
        if max < hist.GetMaximum(): max = hist.GetMaximum();
        hist.SetLineWidth(2)
    
    leg = ROOT.TLegend(0.7,0.7,0.93,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    for i in range(len(names)):
        hists[i].SetLineColor(colors[i]);
        hists[i].SetLineWidth(widths[i]);
        hists[i].SetLineStyle(styles[i]);
        if i == 0 and wGen: leg.AddEntry(hists[i], names[i], "f");
        else: leg.AddEntry(hists[i], names[i], "l");

    if wGen:
        hists[0].SetFillStyle(1001);
        hists[0].SetFillColor(18);
        #hists[0].SetMarkerColor(colors[0]);
        hists[0].SetMarkerStyle(25);
        hists[0].SetLineWidth(0);

    banner = ROOT.TLatex(0.20,0.85,("Anti-kT (R=0.7)"));
    banner.SetNDC()
    banner.SetTextSize(0.035)
    banner2 = ROOT.TLatex(0.20,0.80,("n_{PU} = "+str(options.nPU)));
    banner2.SetNDC()
    banner2.SetTextSize(0.035)
    
    can = ROOT.TCanvas("can"+canname,"can"+canname,1000,800);
    hists[0].SetMaximum( 1.3*max );
    hists[0].SetMinimum( 0 );    
    hists[0].Draw("hist");
    for i in range(1,len(hists)):
        hists[i].Draw("sames");
    leg.Draw();
    banner.Draw();
    banner2.Draw();
    if isLog:
        ROOT.gPad.SetLogy();
        hists[0].SetMinimum( 1 );
    can.SaveAs(directory+"/"+canname+".eps");
    can.SaveAs(directory+"/"+canname+".png");
    can.SaveAs(directory+"/"+canname+".pdf");

    for hist in hists:
        hist.Scale(1./hist.Integral());

def getmax( hs ):
    max = -999;
    for i in range(len(hs)):
        if hs[i] .GetMaximum() > max: max = hs[i].GetMaximum();
    return max;

if __name__ == '__main__':

    #file = ROOT.TFile("../output/outtree_"+str(options.nPU)+".root");
    file = ROOT.TFile("../output/outtree_ttb.root");
    tree_gen = file.Get("tree_gen");
    tree_pf = file.Get("tree_pf");
    tree_pfchs = file.Get("tree_pfchs");
    #tree_pf_tr = file.Get("tree_pf_tr");    
    #tree_pf_cl = file.Get("tree_pf_cl");    
    tree_pf_puppi = file.Get("tree_pf_puppi");    
    
    #types = ["gen","pf","pfchs","pf_tr","pf_cl","pf_puppi"];
    #trees = [tree_gen, tree_pf, tree_pfchs, tree_pf_tr, tree_pf_cl,tree_pf_puppi];
    types = ["GEN","PUPPI","PFlow","PFlowCHS"];
    trees = [tree_gen,tree_pf_puppi,tree_pf,tree_pfchs];
    colors = [1,2,4,6,7,2];
    linestyles = [2,1,1,1,1,2];
    nTypes = len(types);
    
    print "tree_gen.GetEntries() = ",tree_gen.GetEntries();
    print "tree_pf.GetEntries() = ",tree_pf.GetEntries();
    #print "tree_pfchs.GetEntries() = ",tree_pfchs.GetEntries();
    #print "tree_pf_tr.GetEntries() = ",tree_pf_tr.GetEntries();
    #print "tree_pf_cl.GetEntries() = ",tree_pf_cl.GetEntries();
    print "tree_pf_puppi.GetEntries() = ",tree_pf_puppi.GetEntries();

    ## -----------------------------
    ## declare histograms
    h_nJets = [];
    h_jet_mass = [];
    h_jet_mass_trimmed = [];        
    h_jet_pt = [];
    h_jets_pt = [];
    h_jets_eta = [];
    h_jets_response = [];   
    h_jets_response_corr = [];   
    h_jets_mresponse = [];   

    h_sumEt = [];
    h_etMissX = [];
    h_etMissY = [];
    
    h_mass_corr = [];
    for i in range(nTypes):
        
        h_nJets.append(ROOT.TH1F("h_nJets_"+types[i],"; n jets ;count",20,0,20));
        h_jet_mass.append(ROOT.TH1F("h_jet_mass_"+types[i],"; mass (GeV) ;count",20,0,200));
        h_jet_mass_trimmed.append(ROOT.TH1F("h_jet_mass_trimmed_"+types[i],"; trimmed mass (GeV) ;count",20,0,200));
        h_jet_pt.append(ROOT.TH1F("h_jet_pt_"+types[i],"; pT ;count",20,0,800));         

        h_sumEt.append(ROOT.TH1F("h_sumEt_"+types[i],"; #Delta #Sigma Et (GeV) ;count",50,-500,2500));
        h_etMissX.append(ROOT.TH1F("h_etMissX_"+types[i],"; #Delta E_{X}^{miss} (GeV) ;count",50,-200,200));
        h_etMissY.append(ROOT.TH1F("h_etMissY_"+types[i],"; #Delta E_{Y}^{miss} (GeV) ;count",50,-200,200));
        
        h_jets_pt.append(ROOT.TH1F("h_jets_pt"+types[i],"; pT (GeV); count", 100, 20, 720));
        h_jets_eta.append(ROOT.TH1F("h_jets_eta"+types[i],"; eta; count", 30, -5, 5));
        h_jets_response.append(ROOT.TH1F("h_jets_response"+types[i],"; pT - pT_{LV} (GeV); count", 40, -80, 80));
        h_jets_mresponse.append(ROOT.TH1F("h_jets_mresponse"+types[i],"; m - m_{LV} (GeV); count", 40, -80, 80));
        h_jets_response_corr.append(ROOT.TH1F("h_jets_response_corr"+types[i],"; pT/pT_{Gen}; count", 40, -80, 80));

        h_mass_corr.append(ROOT.TH2F("h_mass_corr_"+types[i],"; true mass; reconstructed mass",50,0,200,50,0,200));
                                    
    ## -----------------------------
    ## fill histograms    
    iAreaSubtractors = 1;
    for j in range(nTypes):
        
        print "type = ", j;
        
        for i in range(trees[j].GetEntries()):
            tree_gen.GetEntry(i);
            trees[j].GetEntry(i);

            if i%1000 == 0: print "i = ",i
            #if i > 1000: break;

            curNjets = 0;

            h_sumEt[j].Fill(trees[j].sumEt - tree_gen.sumEt);
            h_etMissX[j].Fill(trees[j].etMissX - tree_gen.etMissX);
            h_etMissY[j].Fill(trees[j].etMissY - tree_gen.etMissY);
            
            for kk in range(trees[j].njets):

                cur4vec = ROOT.TLorentzVector();
                if j > iAreaSubtractors: 
                    cur4vec.SetPtEtaPhiM(trees[j].v_jet_pt_4Vcorr[kk], 
                                        trees[j].v_jet_eta_4Vcorr[kk], 
                                        trees[j].v_jet_phi_4Vcorr[kk], 
                                        trees[j].v_jet_m_4Vcorr[kk] );
                else:
                    cur4vec.SetPtEtaPhiM(trees[j].v_jet_pt[kk], trees[j].v_jet_eta[kk], trees[j].v_jet_phi[kk], trees[j].v_jet_m[kk]);

                if kk == 0:
                    h_jet_mass[j].Fill( cur4vec.M() );
                    h_jet_mass_trimmed[j].Fill( trees[j].v_jet_m_trimmed[0] );            
                    h_jet_pt[j].Fill( cur4vec.Pt() ); 
                    
                if cur4vec.Pt() > 25: 
                    curNjets += 1;
                    h_jets_eta[j].Fill( cur4vec.Eta() );
                    h_jets_pt[j].Fill( cur4vec.Pt() );
                
                # do gen matching
                cur4V = ROOT.TLorentzVector();
                cur4V.SetPtEtaPhiM( trees[j].v_jet_pt[kk], 
                                    trees[j].v_jet_eta[kk], 
                                    trees[j].v_jet_phi[kk], 
                                    trees[j].v_jet_m[kk] );
                cur4Vcorr = ROOT.TLorentzVector();
                cur4Vcorr.SetPtEtaPhiM( trees[j].v_jet_pt_4Vcorr[kk], 
                                    trees[j].v_jet_eta_4Vcorr[kk], 
                                    trees[j].v_jet_phi_4Vcorr[kk], 
                                    trees[j].v_jet_m_4Vcorr[kk] );
                imatch = -1;
                mindr = 0.4;
                for ig in range(tree_gen.njets):
                    curGen4V = ROOT.TLorentzVector();
                    curGen4V.SetPtEtaPhiM(tree_gen.v_jet_pt[ig],tree_gen.v_jet_eta[ig],tree_gen.v_jet_phi[ig],tree_gen.v_jet_m[ig]); 
                    curDR = cur4V.DeltaR(curGen4V);
                    if curDR < mindr: 
                        mindr = cur4V.DeltaR(curGen4V);
                        imatch = ig;
                if imatch >= 0:
                    response = trees[j].v_jet_pt[kk] - tree_gen.v_jet_pt[imatch];
                    response_corr = trees[j].v_jet_pt_4Vcorr[kk] - tree_gen.v_jet_pt[imatch];
                    mresponse = cur4Vcorr.M() - tree_gen.v_jet_m[imatch];
                    #print trees[j].v_jet_pt[kk],tree_gen.v_jet_pt[imatch],response
                    h_jets_response[j].Fill( response );
                    h_jets_response_corr[j].Fill( response_corr );
                    h_jets_mresponse[j].Fill( mresponse );
                                        
            h_nJets[j].Fill(curNjets);

            
    for j in range(trees[0].GetEntries()):
        
        for i in range(nTypes): trees[i].GetEntry(j);
        for i in range(nTypes):
            h_mass_corr[i].Fill(tree_gen.v_jet_m[0],trees[i].v_jet_m[0]);

    ## -----------------------------
    ## Plotting time
    Norm = True;
    makeCanvas(h_nJets,types,"njets",False,Norm,True);
    makeCanvas(h_jet_mass,types,"m_raw",False,Norm,True);
    makeCanvas(h_jet_mass_trimmed,types,"m_trimmed",False,Norm,True);
    makeCanvas(h_jet_pt,types,"pt_corr",False,Norm,True);
    makeCanvas(h_jets_eta,types,"jeteta_corr",False,False,True);
    makeCanvas(h_jets_pt,types,"pts",False,Norm,True);

    massLeg = ["PUPPI","PFlow","PFlowCHS"];
    del h_jets_mresponse[0];
    makeCanvas(h_jets_mresponse,massLeg,"massResponse",False,Norm);
    
    metLeg = ["PUPPI","PFlow","PFlowCHS"];
    del h_sumEt[0];
    del h_etMissX[0];
    del h_etMissY[0];
    makeCanvas(h_sumEt,metLeg,"deltaSumEt",False,Norm);
    makeCanvas(h_etMissX,metLeg,"deltaEXmiss",False,Norm);
    makeCanvas(h_etMissY,metLeg,"deltaEYmiss",False,Norm);

#    makeCanvas(h_jets_response,types,"responses",False,Norm);
#    makeCanvas(h_jets_response_corr,types,"responsesCorr",False,Norm);

    print len(h_jets_response)
    h_jets_response.append( h_jets_response_corr[2] );
    h_jets_response.append( h_jets_response_corr[3] );
    print len(h_jets_response)
    del h_jets_response[0]
    typesClone = ["PUPPI","PFlow (raw)","PFlowCHS (raw)","PFlow (corr.)","PFlowCHS (corr.)"];    
    makeCanvas(h_jets_response,typesClone,"ptResponse",False,Norm);

#    ## -----------------------------
#    ## properties of histograms            
#    for j in range (nTypes):
#        
#        h_nJets[j].Scale(1./h_nJets[j].Integral());
#        h_jet_mass[j].Scale(1./h_jet_mass[j].Integral());
#        h_jet_pt[j].Scale(1./h_jet_pt[j].Integral());
#
#        h_nJets[j].SetLineColor(colors[j]);
#        h_jet_mass[j].SetLineColor(colors[j]);
#        h_jet_pt[j].SetLineColor(colors[j]);
#        h_nJets[j].SetLineWidth(2);
#        h_jet_mass[j].SetLineWidth(2);
#        h_jet_pt[j].SetLineWidth(2);
#        
#        h_nJets[j].SetLineStyle(linestyles[j]);
#        h_jet_mass[j].SetLineStyle(linestyles[j]);
#        h_jet_pt[j].SetLineStyle(linestyles[j]);
#        
#    leg = ROOT.TLegend(0.7,0.7,0.9,0.9);
#    leg.SetBorderSize(0);
#    leg.SetFillStyle(0);
#    for j in range(nTypes):
#        leg.AddEntry(h_nJets[j],types[j],"l");
#
#    # plot
#    can_nJets = ROOT.TCanvas("can_nJets","can_nJets",800,800);
#    h_nJets[0].SetMaximum( 1.2*getmax(h_nJets) );
#    h_nJets[0].Draw();
#    for j in range (1,nTypes): h_nJets[j].Draw("sames");
#    leg.Draw();
#    can_nJets.SaveAs("figs/can_nJets.png");
#
#    can_jet_mass = ROOT.TCanvas("can_jet_mass","can_jet_mass",800,800);
#    h_jet_mass[0].SetMaximum( 1.2*getmax(h_jet_mass) );
#    h_jet_mass[0].Draw();
#    for j in range (1,nTypes): h_jet_mass[j].Draw("sames");
#    leg.Draw();    
#    can_jet_mass.SaveAs("figs/can_jet_mass.png");
#
#    can_jet_pt = ROOT.TCanvas("can_jet_pt","can_jet_pt",800,800);
#    h_jet_pt[0].SetMaximum( 1.2*getmax(h_jet_pt) );
#    h_jet_pt[0].Draw();
#    for j in range (1,nTypes): h_jet_pt[j].Draw("sames");
#    leg.Draw();    
#    can_jet_pt.SaveAs("figs/can_jet_pt.png");

    cantmp = ROOT.TCanvas("cancorr"+str(i),"cancorr"+str(i),800,800);
    drawoption = ["","box","box same","box same"]
    leg = ROOT.TLegend(0.6,0.2,0.95,0.4);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);

    banner = ROOT.TLatex(0.20,0.85,("Anti-kT (R=0.7)"));
    banner.SetNDC()
    banner.SetTextSize(0.035)
    banner2 = ROOT.TLatex(0.20,0.80,("n_{PU} = "+str(options.nPU)));
    banner2.SetNDC()
    banner2.SetTextSize(0.035)

    for i in range(1,nTypes):
        print "corr factor ",i," = ", h_mass_corr[i].GetCorrelationFactor();
        h_mass_corr[i].SetLineColor(colors[i-1]);
        h_mass_corr[i].Draw(drawoption[i]);
        leg.AddEntry(h_mass_corr[i],types[i]+", r = "+str(round(h_mass_corr[i].GetCorrelationFactor(),2)),"l");
    
    leg.Draw();
    banner.Draw();
    banner2.Draw();
    cantmp.SaveAs("plots/mcorr.png");
    cantmp.SaveAs("plots/mcorr.pdf");

        
        