#include "/Users/Phil/crap/tmp/NYStyle.h"
TTree *load(std::string iName) { 
  TFile *lFile = new TFile(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("tree");
  return lTree;
}
void plotProfile(std::string iProfVar,std::string iWeight,std::string iCut,TTree *iTree,double iMax) { 
  TCanvas *lC3 = new TCanvas((iProfVar+"D").c_str(),(iProfVar+"D").c_str(),800,600); 
  TH1F *lH  = new TH1F((iProfVar+"Ring").c_str(),(iProfVar+"Ring").c_str()  ,20,0,iMax); 
  TH1F *lN  = new TH1F((iProfVar+"RNor").c_str(),(iProfVar+"RNor").c_str()  ,20,0,iMax); 
  TH1F *lHP = new TH1F((iProfVar+"RingP").c_str(),(iProfVar+"RingP").c_str(),20,0,iMax); 
  TH1F *lNP = new TH1F((iProfVar+"RNorP").c_str(),(iProfVar+"RNorP").c_str(),20,0,iMax); 
  iTree->Draw((iProfVar+">>"+iProfVar+"Ring").c_str() ,(iWeight+"*(pu == 0)").c_str());
  iTree->Draw((iProfVar+">>"+iProfVar+"RNor").c_str() ,(iCut   +"*(pu == 0)").c_str());
  iTree->Draw((iProfVar+">>"+iProfVar+"RingP").c_str(),(iWeight+"*(pu == 1)").c_str());
  iTree->Draw((iProfVar+">>"+iProfVar+"RNorP").c_str(),(iCut   +"*(pu == 1)").c_str());
  lH ->Divide(lN);  lH ->SetLineColor(kBlue);
  lHP->Divide(lNP); lHP->SetLineColor(kRed); 
  lH ->Draw();
  lHP->Draw("sames");
}
void plot2D(std::string iVar="puppi*ptodRSO",std::string iCut = "(pt > 0.5 )") { 
  Prep();
  TTree *lTree = load("output/OutputTmp.root");
  float lMin  = 1e9;
  float lMax  = -1e9;
  /*
  float lVar = 0; 
  lTree->SetBranchAddress(iVar.c_str(),&lVar);
  for(int i0 = 0; i0 < lTree->GetEntries(); i0++) {
    lTree->GetEntry(i0);
    if(lVar < lMin) lMin = lVar;
    if(lVar > lMax) lMax = lVar;
  }
  */
  //lMax *= 0.5;
  TCanvas *lC0 = new TCanvas("A","A",800,600);
  TH2F* lH2D = new TH2F("Hist","Hist",10,-1.0,1.0,10,-1.0,1.0); 
  TH2F* lN2D = new TH2F("Norm","Norm",10,-1.0,1.0,10,-1.0,1.0); 
  std::string lWeight = iVar + "*" + iCut;
  lTree->Draw("dPhi:dEta>>Hist",(lWeight+"*(pu == 0)").c_str(),"colz");
  lTree->Draw("dPhi:dEta>>Norm",(iCut   +"*(pu == 0)").c_str(),"colz");
  lH2D->Divide(lN2D);
  for(int i0   = 0; i0 < lH2D->GetNbinsX()+1; i0++) { 
    for(int i1 = 0; i1 < lH2D->GetNbinsX()+1; i1++) { 
      if(lMax < lH2D->GetBinContent(i0,i1)) lMax =  lH2D->GetBinContent(i0,i1);
      if(lMin > lH2D->GetBinContent(i0,i1)) lMin =  lH2D->GetBinContent(i0,i1);
    }
  }

  lH2D->GetZaxis()->SetRangeUser(lMin,lMax);
  lH2D->Draw("colz");
  // lC0->SetLogz();
 
  TCanvas *lC1 = new TCanvas("B","B",800,600);
  TH2F* lH2DP = new TH2F("PHist","PHist",10,-1.0,1.0,10,-1.0,1.0); 
  TH2F* lN2DP = new TH2F("PNorm","PNorm",10,-1.0,1.0,10,-1.0,1.0);  
  lTree->Draw("dPhi:dEta>>PHist",(lWeight+"*(pu == 1)").c_str(),"colz");
  lTree->Draw("dPhi:dEta>>PNorm",(iCut   +"*(pu == 1)").c_str(),"colz");  
  lH2DP->Divide(lN2DP);
  lH2DP->GetZaxis()->SetRangeUser(lMin,lMax);
  lH2DP->Draw("colz");
  //lC1->SetLogz();

  TCanvas *lC2 = new TCanvas("C","C",800,600); 
  int lNRings = 3;
  TH1F **pH = new TH1F*[2*lNRings];
  for(int i0 = 0; i0  < lNRings; i0++) { 
    std::stringstream pSH0,pSH1,pSV0,pSV1,pCut0,pCut1; 
    pSH0 << "RingHist" << i0;
    pSH1 << "RingHistPU" << i0;
    pSV0 << iVar << ">>RingHist"   << i0;   
    pSV1 << iVar << ">>RingHistPU" << i0;
    pCut0 << "(" << i0*0.15 << " < dr && dr < " << i0*(0.15) + 0.15 << " )*" << iCut;    
    pCut1 << "(" << i0*0.15 << " < dr && dr < " << i0*(0.15) + 0.15 << " )*" << iCut;
    pH[2.*i0]   = new TH1F(pSH0.str().c_str(),pSH0.str().c_str(),30,lMin,lMax*5.);
    pH[2.*i0+1] = new TH1F(pSH1.str().c_str(),pSH1.str().c_str(),30,lMin,lMax*5.);
    lTree->Draw(pSV0.str().c_str(),(pCut0.str()+"*(pu == 0)").c_str()); 
    lTree->Draw(pSV1.str().c_str(),(pCut1.str()+"*(pu == 1)").c_str());
    pH[2.*i0]  ->Scale(1./pH[2.*i0]  ->Integral());
    pH[2.*i0+1]->Scale(1./pH[2.*i0+1]->Integral());
    pH[2.*i0]  ->SetLineColor(kBlue+i0); 
    pH[2.*i0+1]->SetLineColor(kRed +i0); 
    pH[2.*i0]  ->SetLineStyle(1+i0);
    pH[2.*i0+1]->SetLineStyle(1+i0);
  }
  pH[0]->Draw();
  for(int i0 = 1; i0 < lNRings; i0++) pH[i0]->Draw("sames");  
  plotProfile("dr",lWeight,iCut,lTree,0.8);
  plotProfile("pt",lWeight,iCut,lTree,5.);
}
