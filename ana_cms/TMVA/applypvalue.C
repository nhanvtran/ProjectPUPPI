#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <iostream>

float fMVA    = 0; 
int   fPFType = 0;
float fGenPt  = 0; 
float fPt     = 0; 
int isGen() { 
  double pVal = fGenPt/fPt;
  if(pVal < 0)               return -1;
  if(pVal > 0 && pVal < 0.5) return  0;
  return 1;
}
void meanRMS(float &iMean,float &iRMS,std::vector<float> &iVMean) { 
  iMean = 0; iRMS = 0;
  for(unsigned int i0 = 0; i0 < iVMean.size(); i0++) { 
    iMean += iVMean[i0];
  } 
  iMean/=iVMean.size();
  for(unsigned int i0 = 0; i0 < iVMean.size(); i0++) { 
    iRMS += (iVMean[i0]-iMean)*(iVMean[i0]-iMean);
  } 
  iRMS=sqrt(iRMS/iVMean.size());
}
void computeMeanRMS(float *iMean,float *iRMS,TTree *iTree) { 
  std::vector<float> lVals0;
  std::vector<float> lVals1;
  std::vector<float> lVals2;
  std::vector<float> lVals0PU;
  std::vector<float> lVals1PU;
  std::vector<float> lVals2PU;
  for(int i0 = 0; i0 < iTree->GetEntries(); i0++) { 
    iTree->GetEntry(i0);
    int lGen = isGen();
    if(lGen == 0) continue;
    if(lGen == 1) {
      if(fPFType > 0 && fPFType < 4) lVals0.push_back(fMVA);
      if(fPFType == 4)               lVals1.push_back(fMVA);
      if(fPFType == 5)               lVals2.push_back(fMVA);
      continue;
    }
    if(fPFType > 0 && fPFType < 4)   lVals0PU.push_back(fMVA);
    if(fPFType == 4)                 lVals1PU.push_back(fMVA);
    if(fPFType == 5)                 lVals2PU.push_back(fMVA);
  }
  iMean[0] = 0; iRMS[0] = 0; meanRMS(iMean[0],iRMS[0],lVals0);
  iMean[1] = 0; iRMS[1] = 0; meanRMS(iMean[1],iRMS[1],lVals0PU);
  iMean[2] = 0; iRMS[2] = 0; meanRMS(iMean[2],iRMS[2],lVals1);
  iMean[3] = 0; iRMS[3] = 0; meanRMS(iMean[3],iRMS[3],lVals1PU);
  iMean[4] = 0; iRMS[4] = 0; meanRMS(iMean[4],iRMS[4],lVals2);
  iMean[5] = 0; iRMS[5] = 0; meanRMS(iMean[5],iRMS[5],lVals2PU);
  for(int i0 = 0; i0 < 6; i0++) { 
    std::cout << " Mean : " << i0 << " -- " << iMean[i0] << " - " << iRMS[i0] << std::endl;
  }
}
float chi2(int iId,float *iMean,float *iRMS) { 
  int lId = 0; 
  if(fPFType == 4)  lId = 2;
  if(fPFType == 5)  lId = 4;
  lId += iId;
  float lMean = iMean[lId];
  float lRMS  = iRMS [lId];
  float lChi2 = (fMVA-lMean)/lRMS;
  return lChi2*lChi2;
}
void applypvalue(std::string iName="BDTOutput.root")  { 
  TFile *lFile = new TFile(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Result"); 
  fMVA    = 0; lTree->SetBranchAddress("bdt"   ,&fMVA);
  fPFType = 0; lTree->SetBranchAddress("pftype",&fPFType);
  fGenPt  = 0; lTree->SetBranchAddress("genpt" ,&fGenPt);
  fPt     = 0; lTree->SetBranchAddress("pt"    ,&fPt);

  float *lMean = new float[6];
  float *lRMS  = new float[6];
  computeMeanRMS(lMean,lRMS,lTree);
  TFile *lOFile = new TFile("Output.root","RECREATE");
  TTree *lOTree  = (TTree*) lTree->CloneTree(0);
  float  lChi2PV = 0; lOTree->Branch("chi2pv",&lChi2PV,"lChi2PV/F");
  float  lChi2PU = 0; lOTree->Branch("chi2pu",&lChi2PU,"lChi2PU/F");
  float  lProb   = 0; lOTree->Branch("prob"  ,&lProb  ,"lProb/F"  );
  for(int i0 = 0; i0 < lTree->GetEntries(); i0++) { 
    lTree->GetEntry(i0);
    lChi2PV = chi2(0,lMean,lRMS);
    lChi2PU = chi2(1,lMean,lRMS);
    float pP1 = TMath::Prob(lChi2PV,1.);
    float pP2 = TMath::Prob(lChi2PU,1.);
    lProb = pP1*(1-pP2);
    lProb *= 2.; if(lProb > 1) lProb = 1;
    //if(fPFType == 1) std::cout << " --> " << lChi2PV << " -- " << pP1 << " -- " << lChi2PU << " -- " << pP2 << " -- " << -2*log(pP1/pP2) << " -- " << lProb << std::endl;
    if(fPFType > 5) {lProb = 1;}
    lOTree->Fill();
  }
  lOTree->Write();
}
