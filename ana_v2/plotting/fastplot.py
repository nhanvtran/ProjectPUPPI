#! /usr/bin/env python
import os
import glob
import math
import array
import sys
import time

from array import array

import ROOT
from ROOT import gROOT, gStyle, gSystem, TLatex

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadTopMargin(0.09);
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPalette(1);
ROOT.gStyle.SetErrorX(0.5);

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--type',action="store",type="string",dest="type",default="gg")
parser.add_option('--ptbin',action="store",type="int",dest="ptbin",default=500)

(options, args) = parser.parse_args()


############################################################
def makeROCFromHisto(hsig,hbkg,LtoR):
    
    nbins = hsig.GetNbinsX();
    binsize = hsig.GetBinWidth(1);
    lowedge = hsig.GetBinLowEdge(1);
    
    #print "lowedge: ",lowedge
    
    hsigIntegral = hsig.Integral();
    hbkgIntegral = hbkg.Integral();
    
    xval = array('d', [])
    yval = array('d', [])
    ctr = 0;
    for i in range(1,nbins+1):
        
        effBkg = 0;
        effSig = 0;
        
        if LtoR: effBkg = hbkg.Integral( i, nbins )/hbkgIntegral;
        else: effBkg = hbkg.Integral( 1, i )/hbkgIntegral;
        
        if LtoR: effSig = hsig.Integral( i, nbins )/hsigIntegral;
        else: effSig = hsig.Integral( 1, i )/hsigIntegral;
        
        #print "cut: ",(lowedge+(i-1)*binsize),"effBkg: ", effBkg, ", effSig: ", effSig;
        
        xval.append( effSig );
        yval.append( 1-effBkg );
        ctr = ctr + 1;
    
    #print nbins, "and ", ctr
    tg = ROOT.TGraph( nbins, xval, yval );
    return tg;    
    
def makeCanvas(hists, names, canname):
    
    directory = "plots";
    colors = [2,4,1,6,7];
    
    max = -999.;
    for hist in hists:
        if max < hist.GetMaximum(): max = hist.GetMaximum();

    leg = ROOT.TLegend(0.7,0.7,0.9,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    for i in range(len(names)):
        hists[i].SetLineColor(colors[i]);
        leg.AddEntry(hists[i], names[i], "l");

    can = ROOT.TCanvas("can"+canname,"can"+canname,1000,800);
    hists[0].SetMaximum( 1.2*max );
    hists[0].SetMinimum( 0 );    
    hists[0].Draw();    
    for i in range(1,len(hists)):
        hists[i].Draw("sames");
    leg.Draw();
    can.SaveAs(directory+"/"+canname+".eps");
    can.SaveAs(directory+"/"+canname+".png");

    for hist in hists:
        hist.Scale(1./hist.Integral());
        
def makeSingleHistCanvas( h_2D ):
    
    directory = "figs_bin"+str(bin);

    cant = ROOT.TCanvas("cant","cant",800,800);
    h_2D.Draw("colz");
    cant.SaveAs(directory+"/"+h_2D.GetName()+".png");
    cant.SaveAs(directory+"/"+h_2D.GetName()+".eps");
    
def drawROCs(rocs, legtext, name):
        
    legTextSize=0.037;
        
    colors = [1,2,4,6,7,3,8,10]

    banner = TLatex(0.25,0.92,("Daniele Preliminary Simulation, #sqrt{s} = 8 TeV, Z+jets"));
    banner.SetNDC()
    banner.SetTextSize(0.035)

    leg = ROOT.TLegend( 0.2, 0.20, 0.6, 0.45 );
    leg.SetBorderSize( 0 );
    leg.SetFillStyle( 0 );
    leg.SetNColumns(1);
    for i in range(len(rocs)):
        leg.AddEntry( rocs[i], legtext[i], "l" );
    
    canRoc = ROOT.TCanvas("canRoc","canRoc",800,800);    
    canRoc.cd();
    hrl = canRoc.DrawFrame(0,0,1.0,1.0);
    hrl.GetXaxis().SetTitle("#epsilon_{sig}");
    hrl.GetYaxis().SetTitle("1 - #epsilon_{bkg}");    
    for i in range(len(rocs)):
        rocs[i].SetLineColor( colors[i] );
        rocs[i].Draw();
    leg.Draw();
    banner.Draw();
        
    canRoc.SaveAs( "plots/"+name+".eps" );
    canRoc.SaveAs( "plots/"+name+".png" );
    canRoc.SaveAs( "plots/"+name+".pdf" );    
    del canRoc;
    
############################################################

if __name__ == '__main__':

    file = ROOT.TFile("../output/outtree_80.root");
    tree = file.Get("tree_particles");
    
    hlo = -1;
    hhi =  1;
    
    weights_PU = ROOT.TH1F("weights_PU","; weight; count", 100, hlo, hhi);
    weights_LV = ROOT.TH1F("weights_LV","; weight; count", 100, hlo, hhi);        

    weights_PU_chLV = ROOT.TH1F("weights_PU_chLV","; weight; count", 100, hlo, hhi);
    weights_LV_chLV = ROOT.TH1F("weights_LV_chLV","; weight; count", 100, hlo, hhi);        

    weights_PU_chPU = ROOT.TH1F("weights_PU_chPU","; weight; count", 100, hlo, hhi);
    weights_LV_chPU = ROOT.TH1F("weights_LV_chPU","; weight; count", 100, hlo, hhi);        

    weights_PU_ch = ROOT.TH1F("weights_PU_ch","; weight; count", 100, 2*hlo, 2*hhi);
    weights_LV_ch = ROOT.TH1F("weights_LV_ch","; weight; count", 100, 2*hlo, 2*hhi);        

    weights_PU_chLV_shape = ROOT.TH1F("weights_PU_chLV_shape","; weight; count", 100, 2*hlo, 2*hhi);
    weights_LV_chLV_shape = ROOT.TH1F("weights_LV_chLV_shape","; weight; count", 100, 2*hlo, 2*hhi);        

    weights_PU_all = ROOT.TH1F("weights_PU_all","; weight; count", 100, 2*hlo, 2*hhi);
    weights_LV_all = ROOT.TH1F("weights_LV_all","; weight; count", 100, 2*hlo, 2*hhi);        
    
    for i in range(tree.GetEntries()):
        tree.GetEntry(i);
        cur_pt = math.sqrt(tree.p_px*tree.p_px + tree.p_py*tree.p_py);
        if cur_pt > 1:
            if tree.p_isPU == 0: 
                weights_LV.Fill( tree.p_puppiW );
                weights_LV_chLV.Fill( tree.p_puppiW_chLV );                
                weights_LV_chPU.Fill( tree.p_puppiW_chPU ); 
                weights_LV_ch.Fill( tree.p_puppiW_chLV - tree.p_puppiW_chPU );                                
                weights_LV_chLV_shape.Fill( tree.p_puppiW_chLV + tree.p_puppiW );                
                weights_LV_all.Fill( tree.p_puppiW + tree.p_puppiW_chLV - tree.p_puppiW_chPU );                                
            if tree.p_isPU == 1: 
                weights_PU.Fill( tree.p_puppiW );
                weights_PU_chLV.Fill( tree.p_puppiW_chLV );                
                weights_PU_chPU.Fill( tree.p_puppiW_chPU );                                
                weights_PU_ch.Fill( tree.p_puppiW_chLV - tree.p_puppiW_chPU );                                
                weights_PU_chLV_shape.Fill( tree.p_puppiW_chLV + tree.p_puppiW );                
                weights_PU_all.Fill( tree.p_puppiW + tree.p_puppiW_chLV - tree.p_puppiW_chPU );                                
                                 
    makeCanvas( [weights_PU,weights_LV], ["PU shape","LV shape"], "puppi_shape" );
    makeCanvas( [weights_PU_chLV,weights_LV_chLV], ["PU chLV","LV chLV"], "puppi_chLV" );
    makeCanvas( [weights_PU_chPU,weights_LV_chPU], ["PU chPU","LV chPU"], "puppi_chPU" );
    makeCanvas( [weights_PU_ch,weights_LV_ch], ["PU ch","LV ch"], "puppi_ch" );
    makeCanvas( [weights_PU_chLV_shape,weights_LV_chLV_shape], ["PU chLV+shape","LV chLV+shape"], "puppi_chLV_shape" );
    makeCanvas( [weights_PU_all,weights_LV_all], ["PU all","LV all"], "puppi_all" );        

    roc_shape = makeROCFromHisto(weights_PU,weights_LV, True);
    roc_chLV = makeROCFromHisto(weights_PU_chLV,weights_LV_chLV, True);
    roc_chPU = makeROCFromHisto(weights_PU_chPU,weights_LV_chPU, True);
    roc_ch = makeROCFromHisto(weights_PU_ch,weights_LV_ch, True);
    roc_chLV_shape = makeROCFromHisto(weights_PU_chLV_shape,weights_LV_chLV_shape, True);
    roc_all = makeROCFromHisto(weights_PU_all,weights_LV_all, True);
    
    rocs = [roc_shape,roc_chLV,roc_chPU,roc_ch,roc_chLV_shape,roc_all];
    rocNames = ["shape","chLV","chPU","ch","chLV+shape","shape+chLV-chPU"];
    
    drawROCs(rocs,rocNames,"rocs");
                                                                      
#    can1 = ROOT.TCanvas("can1","can1",1000,800);
#    weights_PU.SetMinimum(1);
#    weights_PU.Draw();
#    weights_LV.Draw("sames");   
#    #ROOT.gPad.SetLogy();
#    can1.SaveAs("test_shape.png");     
#    
#    can2 = ROOT.TCanvas("can2","can2",1000,800);
#    weights_PU_chLV.SetMinimum(1);
#    weights_PU_chLV.Draw();
#    weights_LV_chLV.Draw("sames");   
#    #ROOT.gPad.SetLogy();
#    can2.SaveAs("test_chLV.png");     
#
#    can3 = ROOT.TCanvas("can3","can3",1000,800);
#    weights_PU_chPU.SetMinimum(1);
#    weights_PU_chPU.Draw();
#    weights_LV_chPU.Draw("sames");   
#    #ROOT.gPad.SetLogy();
#    can3.SaveAs("test_chPU.png");     
#
#    can4 = ROOT.TCanvas("can4","can4",1000,800);
#    weights_PU_ch.SetMinimum(1);
#    weights_PU_ch.Draw();
#    weights_LV_ch.Draw("sames");   
#    #ROOT.gPad.SetLogy();
#    can4.SaveAs("test_ch.png");     
#
#    can5 = ROOT.TCanvas("can5","can5",1000,800);
#    weights_PU_all.SetMinimum(1);
#    weights_PU_all.Draw();
#    weights_LV_all.Draw("sames");   
#    #ROOT.gPad.SetLogy();
#    can5.SaveAs("test_all.png");     
