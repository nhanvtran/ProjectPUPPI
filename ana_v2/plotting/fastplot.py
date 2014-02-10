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
ROOT.gStyle.SetPadRightMargin(0.16);
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
    
def makeCanvas(hists, names, canname, isLog=False):
    
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
    if isLog: 
        ROOT.gPad.SetLogy();
        hists[0].SetMinimum( 1 );
    can.SaveAs(directory+"/"+canname+".eps");
    can.SaveAs(directory+"/"+canname+".png");


    for hist in hists:
        hist.Scale(1./hist.Integral());
        
def makeSingleHistCanvas( h_2D ):
    
    directory = "plots";

    h_2D.Scale( 1./h_2D.Integral() );
    
    cant = ROOT.TCanvas("cant","cant",1000,800);
    h_2D.Draw("colz");
    ROOT.gPad.SetLogz();
    cant.SaveAs(directory+"/"+h_2D.GetName()+".png");
    cant.SaveAs(directory+"/"+h_2D.GetName()+".eps");
    
def drawROCs(rocs, legtext, name):
        
    legTextSize=0.037;
        
    colors = [1,2,4,6,7,3,5,10]

    banner = TLatex(0.25,0.92,("Daniele Preliminary Simulation, #sqrt{s} = 8 TeV, Z+jets"));
    banner.SetNDC()
    banner.SetTextSize(0.035)

    leg = ROOT.TLegend( 0.2, 0.20, 0.6, 0.45 );
    leg.SetBorderSize( 0 );
    leg.SetFillStyle( 0 );
    leg.SetNColumns(1);
    for i in range(len(rocs)):
        rocs[i].SetLineWidth(2);
        leg.AddEntry( rocs[i], legtext[i], "l" );
    
    canRoc = ROOT.TCanvas("canRoc","canRoc",1000,800);    
    canRoc.cd();
    hrl = canRoc.DrawFrame(0.0,0.0,1.0,1.0);
    hrl.GetXaxis().SetTitle("#epsilon_{sig}");
    hrl.GetYaxis().SetTitle("1 - #epsilon_{bkg}");    
    for i in range(len(rocs)):
        rocs[i].SetLineColor( colors[i] );
        rocs[i].Draw();
    leg.Draw();
    #banner.Draw();
#    ROOT.gPad.SetLogy();
#    ROOT.gPad.SetLogx();
            
    canRoc.SaveAs( "plots/"+name+".eps" );
    canRoc.SaveAs( "plots/"+name+".png" );
    canRoc.SaveAs( "plots/"+name+".pdf" );    
    del canRoc;
    
def findEqualTailProbabilities(h1,h2):

    tailProbDiff = 9999;
    iEqual = -1.;
    nbinsX = h1.GetNbinsX();
    
    nTotal1 = h1.GetEntries();
    nTotal2 = h2.GetEntries();
    
    for i in range(1,nbinsX+1):
        int1 = h1.Integral(1,i);    
        int2 = h2.Integral(i,nbinsX+1);        
        print "val = ", h1.GetBinCenter( i ), "int1 = ", int1, ", int2 = ", int2, ", tailProbDiff = ", tailProbDiff
        if math.fabs(int1-int2) < tailProbDiff: 
            tailProbDiff = math.fabs(int1-int2)
            iEqual = i;
    
    return h1.GetBinCenter( iEqual );

############################################################

if __name__ == '__main__':

    file = ROOT.TFile("../output/outtree_80.root");
    tree = file.Get("tree_particles");
    
    hlo =  -0.1;
    hhi =  1.1;
    nbins = 100;
#    
#    weights_PU_all = ROOT.TH1F("weights_PU_all",";weights;count",100,-10,10);
#    weights_LV_all = ROOT.TH1F("weights_LV_all",";weights;count",100,-10,10);
#
#    weights_PU_chLV = ROOT.TH1F("weights_PU_chLV",";weights;count",100,-100,100);
#    weights_LV_chLV = ROOT.TH1F("weights_LV_chLV",";weights;count",100,-100,100);
#
#        
#    for i in range(tree.GetEntries()):
#        tree.GetEntry(i);
#        cur_pt = math.sqrt(tree.p_px*tree.p_px + tree.p_py*tree.p_py);
#
#        if cur_pt > 1 and tree.p_isPU == 0:
##            print "tree.p_puppiW = ", tree.p_puppiW
#            weights_LV_all.Fill( tree.p_puppiW );
#            weights_LV_chLV.Fill( tree.p_puppiW_chLV );
#
#        if cur_pt > 1 and tree.p_isPU == 1:
#            weights_PU_all.Fill( tree.p_puppiW );
#            weights_PU_chLV.Fill( tree.p_puppiW_chLV );
#                                                                            
#    makeCanvas( [weights_LV_all,weights_PU_all], ["LV #beta1","PU #beta1"], "puppi_shape" );
#    makeCanvas( [weights_LV_chLV,weights_PU_chLV], ["LV #beta2","PU #beta2"], "puppi_chLV" );

    # for ch only
    weights_all_pfchs = ROOT.TH1F("weights_pfchs","; weight; count", nbins, hlo, hhi);
    weights_PU_pfchs = ROOT.TH1F("weights_PU_pfchs","; weight; count", nbins, hlo, hhi);
    weights_LV_pfchs = ROOT.TH1F("weights_LV_pfchs","; weight; count", nbins, hlo, hhi);        
    weights_nePU_pfchs = ROOT.TH1F("weights_nePU_pfchs","; weight; count", nbins, hlo, hhi);
    weights_neLV_pfchs = ROOT.TH1F("weights_neLV_pfchs","; weight; count", nbins, hlo, hhi);        

    weights_all_chLV = ROOT.TH1F("weights_chLV","; weight; count", nbins, hlo, hhi);
    weights_PU_chLV = ROOT.TH1F("weights_PU_chLV","; weight; count", nbins, hlo, hhi);
    weights_LV_chLV = ROOT.TH1F("weights_LV_chLV","; weight; count", nbins, hlo, hhi);        
    weights_nePU_chLV = ROOT.TH1F("weights_nePU_chLV","; weight; count", nbins, hlo, hhi);
    weights_neLV_chLV = ROOT.TH1F("weights_neLV_chLV","; weight; count", nbins, hlo, hhi);        

    weights_all_all = ROOT.TH1F("weights_all_all","; weight; count", nbins, hlo, hhi);
    weights_PU_all = ROOT.TH1F("weights_PU_all","; weight; count", nbins, hlo, hhi);
    weights_LV_all = ROOT.TH1F("weights_LV_all","; weight; count", nbins, hlo, hhi);        
    weights_nePU_all = ROOT.TH1F("weights_nePU_all","; weight; count", nbins, hlo, hhi);
    weights_neLV_all = ROOT.TH1F("weights_neLV_all","; weight; count", nbins, hlo, hhi);        

    weights_LV_allVschLV = ROOT.TH2F("weights_LV_allVschLV",";shape;LV",nbins,hlo,hhi,nbins,hlo,hhi);
    weights_PU_allVschLV = ROOT.TH2F("weights_PU_allVschLV",";shape;LV",nbins,hlo,hhi,nbins,hlo,hhi);
    
    for i in range(tree.GetEntries()):
        tree.GetEntry(i);
        cur_pt = math.sqrt(tree.p_px*tree.p_px + tree.p_py*tree.p_py);
        
        if cur_pt > 1:
            weights_all_pfchs.Fill( tree.p_puppiW_pfchs);
            weights_all_chLV.Fill( tree.p_puppiW_chLV );
            weights_all_all.Fill( tree.p_puppiW_all );            

        if cur_pt > 1 and tree.p_isPU == 0 and tree.p_isCH == 0:
            weights_neLV_pfchs.Fill( tree.p_puppiW_pfchs);
            weights_neLV_chLV.Fill( tree.p_puppiW_chLV );
            weights_neLV_all.Fill( tree.p_puppiW_all );                        

        if tree.p_isCH <= 1 and tree.p_isPU == 0 and cur_pt > 1:
            weights_LV_pfchs.Fill( tree.p_puppiW_pfchs);
            weights_LV_chLV.Fill( tree.p_puppiW_chLV );
            weights_LV_all.Fill( tree.p_puppiW_all );   

            weights_LV_allVschLV.Fill( tree.p_puppiW_all, tree.p_puppiW_chLV );
                                                                              
        if tree.p_isPU == 1 and tree.p_isCH == 0 and cur_pt > 1: 
            weights_nePU_pfchs.Fill( tree.p_puppiW_pfchs);
            weights_nePU_chLV.Fill( tree.p_puppiW_chLV );                    
            weights_nePU_all.Fill( tree.p_puppiW_all );                    
       
        if tree.p_isCH <= 1 and tree.p_isPU == 1 and cur_pt > 1:
            weights_PU_pfchs.Fill( tree.p_puppiW_pfchs);
            weights_PU_chLV.Fill( tree.p_puppiW_chLV );
            weights_PU_all.Fill( tree.p_puppiW_all );

            weights_PU_allVschLV.Fill( tree.p_puppiW_all, tree.p_puppiW_chLV );    

#    cutOn = findEqualTailProbabilities(weights_LV_chLV,weights_PU_chLV);
#    print "cutOn = ", cutOn;

    makeCanvas( [weights_LV_pfchs,weights_PU_pfchs],   ["LV #Sigma pfchs", "PU #Sigma pfchs"], "puppi_pfchs", True );        
    makeCanvas( [weights_LV_chLV,weights_PU_chLV], ["LV #Sigma chLV","PU #Sigma chLV"], "puppi_chLV", True );        
    makeCanvas( [weights_LV_all,weights_PU_all], ["LV #Sigma all","PU #Sigma all"], "puppi_all", True );        

    makeSingleHistCanvas(weights_LV_allVschLV);
    makeSingleHistCanvas(weights_PU_allVschLV);

    roc_pfchs = makeROCFromHisto(weights_LV_pfchs,weights_PU_pfchs, True);
    roc_chLV = makeROCFromHisto(weights_LV_chLV,weights_PU_chLV, True);
    roc_all = makeROCFromHisto(weights_LV_all,weights_PU_all, True);
#    roc_ch = makeROCFromHisto(weights_PU_ch,weights_LV_ch, True);
#    roc_chLV_shape = makeROCFromHisto(weights_PU_chLV_shape,weights_LV_chLV_shape, True);
#    roc_chLV_shape_prod = makeROCFromHisto(weights_PU_chLV_shape_prod,weights_LV_chLV_shape_prod, False);
#    roc_all = makeROCFromHisto(weights_PU_all,weights_LV_all, True);
#    
    rocs = [roc_pfchs,roc_chLV,roc_all];
    rocNames = ["pfchs","chLV","all"];
    #rocs = [roc_chLV]
    #rocNames = ["puppi tracks"];

    
    drawROCs(rocs,rocNames,"rocs");
