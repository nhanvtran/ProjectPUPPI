std::vector<std::string> fMap;
std::vector<std::string> fVar;
void load() { 
fMap.push_back("ptc");     fVar.push_back("#Sigma p_{T}    (GeV)");
fMap.push_back("dR");      fVar.push_back("#Sigma #Delta R (GeV)");
fMap.push_back("ptdR");    fVar.push_back("#Sigma log(p_{T} #Delta R) ");
fMap.push_back("puppi");   fVar.push_back("#Sigma log(p_{T}/sqrt(#Delta R)) ");
 fMap.push_back("ptodR");   fVar.push_back("#Sigma log(p_{T}/#Delta R) (GeV)");
fMap.push_back("ptodRS");  fVar.push_back("#Sigma log(p_{T}/#Delta R/#Sigma) (GeV)");
fMap.push_back("ptodRSO");  fVar.push_back("#Sigma log(p_{T}/#Delta R)/#Sigma (GeV)");

fMap.push_back("ptc_lv");   fVar.push_back("leading #Sigma p_{T}    (GeV)");
fMap.push_back("dR_lv");   fVar.push_back("leading #Sigma #Delta R (GeV)");
fMap.push_back("ptdR_lv");  fVar.push_back("leading #Sigma log(p_{T} #Delta R) ");
fMap.push_back("puppi_lv"); fVar.push_back("leading #Sigma log(p_{T}/sqrt(#Delta R)) ");
 fMap.push_back("ptodR_lv"); fVar.push_back("leading #Sigma log(p_{T}/#Delta R) (GeV)");
fMap.push_back("ptodRSO_lv");fVar.push_back("leading #Sigma log(p_{T}/#Delta R)/#Sigma (GeV)");

fMap.push_back("pt_pu");   fVar.push_back("pu #Sigma p_{T}    (GeV)");
fMap.push_back("dR_pu");   fVar.push_back("pu #Sigma #Delta R (GeV)");
fMap.push_back("ptdR_pu");  fVar.push_back("pu #Sigma log(p_{T} #Delta R) ");
fMap.push_back("puppi_pu"); fVar.push_back("pu #Sigma log(p_{T}/sqrt(#Delta R)) ");
 fMap.push_back("ptodR_pu"); fVar.push_back("pu #Sigma log(p_{T}/#Delta R/#Sum p_{T}) (GeV)");
fMap.push_back("ptodRS_pu"); fVar.push_back("pu #Sigma log(p_{T}/#Delta R/#Sum p_{T}) (GeV)");
fMap.push_back("ptodRSO_pu");fVar.push_back("pu #Sigma log(p_{T}/#Delta R)/#Sigma (GeV)");
}
std::string getName(std::string iVar) { 
  int iId = 0; 
  for(int i0 = 0; i0 < fMap.size(); i0++) { 
    if(iVar == fMap[i0]) iId = i0;
  }
  return fVar[iId];
}

TTree *load(std::string iName) { 
  TFile *lFile = new TFile(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("tree");
  return lTree;
}
TGraph *roc(TH1F* iH0,TH1F *iH1,bool iSames=false) { 
  int lN = iH0->GetNbinsX(); 
  double *lX = new double[lN];
  double *lY = new double[lN];
  bool lFirst = true;
  for(int i0 = 1; i0 < iH0->GetNbinsX()+1; i0++) { 
    lX[i0-1] = (iH0->Integral(lN-i0,1e9))/iH0->Integral(0,1e9); 
    lY[i0-1] = (iH1->Integral(lN-i0,1e9))/iH1->Integral(0,1e9); 
    if(lY[i0-1] > 0.85 && lFirst) {cout << "---> Back at 85% : " << lX[i0-1] << endl; lFirst = false;}
  }
  TGraph *lGraph = new TGraph(lN,lX,lY); lGraph->SetLineColor(kRed); lGraph->SetLineWidth(3);
  if(iSames) lGraph->SetLineColor(kBlue);
  lGraph->SetTitle("");
  lGraph->GetXaxis()->SetTitle("#epsilon_{back}");
  lGraph->GetYaxis()->SetTitle("#epsilon_{sig}");
  if(!iSames) lGraph->Draw("al");
  if(iSames) lGraph->Draw("l");
}
void makeSimplePage(std::string iVar="puppi*ptodRSO",std::string iCut="(pt > 1. && dr > 0.2 )") { 
  TTree *lTree = load("output/OutputTmp.root");
  float lMin  = 1e9;
  float lMax  = -1e9;
  float lVar = 0; 
  /*
  lTree->SetBranchAddress(iVar.c_str(),&lVar);
  for(int i0 = 0; i0 < lTree->GetEntries(); i0++) {
    lTree->GetEntry(i0);
    if(lVar < lMin) lMin = lVar;
    if(lVar > lMax) lMax = lVar;
  }
  */
  lMin = -5; lMax = 200;
  TH1F* lH    = new TH1F("A","A",100,lMin,lMax); lH   ->SetLineColor(kBlue);  lH   ->SetLineWidth(2);
  TH1F* lHPU  = new TH1F("B","B",100,lMin,lMax); lHPU ->SetLineColor(kRed);   lHPU ->SetLineWidth(2);
  TH1F* lH1   = new TH1F("C","C",100,lMin,lMax); lH1  ->SetLineColor(kBlue);  lH1  ->SetLineWidth(2); lH1  ->SetLineStyle(kDashed); 
  TH1F* lHPU1 = new TH1F("D","D",100,lMin,lMax); lHPU1->SetLineColor(kRed);   lHPU1->SetLineWidth(2); lHPU1->SetLineStyle(kDashed);
  lTree->Draw((iVar+">>A").c_str()    ,(iCut+"*(pu == 0)").c_str());
  lTree->Draw((iVar+">>B").c_str()    ,(iCut+"*(pu >  0)").c_str());
  lTree->Draw((iVar+">>C").c_str()    ,(iCut+"*(pu == 0 && charge == 0)").c_str());
  lTree->Draw((iVar+">>D").c_str()    ,(iCut+"*(pu >  0 && charge == 0)").c_str());
  lH->GetXaxis()->SetTitle((getName(iVar)).c_str());
  TCanvas *lCan = new TCanvas("A","A",800,600);
  lH   ->Draw();
  lHPU ->Draw("sames");
  lH1  ->Draw("sames");
  lHPU1->Draw("sames");
  //lCan->SaveAs((iVar+".png").c_str());
  TCanvas *lCan = new TCanvas("B","B",800,600);
  roc(lHPU1,lH1);
}
void makePage() { 
  load();
  makeSimplePage("ptodRSO+puppi");
  //for(int i0 = 0; i0 < fMap.size(); i0++) makeSimplePage(fMap[i0]);
}
