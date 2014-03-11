#include "puppiTMVAContainer.hh"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "fastjet/internal/base.hh"
#include "TH2F.h"
#include "TMath.h"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

puppiTMVAContainer::puppiTMVAContainer(std::vector<RecoObj> &inParticles,std::string iWeight) { 
  TMVA::Tools::Instance();
  fReader      = new TMVA::Reader( "!Color:!Silent" );    
  fPt       = 0; fReader->AddVariable("pt"                 , &fPt);
  fEta      = 0; fReader->AddVariable("eta"                , &fEta);
  fTk       = 0; fReader->AddVariable("tk"                 , &fTk);
  fPVtx     = 0; fReader->AddVariable("vtx"                , &fPVtx);
  fVId      = 0; fReader->AddVariable("vid"                , &fVId);
  fDepth    = 0; fReader->AddVariable("depth"              , &fDepth);
  fTime     = 0; fReader->AddVariable("time"               , &fTime);
  fPFType   = 0; fReader->AddVariable("pfType"             , &fPFType);
  fPuppi    = 0; fReader->AddVariable("ptodR"              , &fPuppi);
  fPuppiO   = 0; fReader->AddVariable("ptodRSO"            , &fPuppiO);
  fPuppiLV  = 0; fReader->AddVariable("ptodR_lv"           , &fPuppiLV);
  fPuppiOLV = 0; fReader->AddVariable("ptodRSO_lv"         , &fPuppiOLV);
  fReader     ->BookMVA("BDT",iWeight.c_str());
  refresh(inParticles); 
}
void puppiTMVAContainer::refresh(std::vector<RecoObj> &inParticles) { 
  _recoParticles.resize(0); 
  _pfParticles.resize(0);
  _genParticles.resize(0);
  _pfchsParticles.resize(0);    
  _chargedPV.resize(0); 
  _chargedNoPV.resize(0);
  _vals.resize(0); 

  _recoParticles = inParticles;
 for (unsigned int i = 0; i < inParticles.size(); i++){
    fastjet::PseudoJet curPseudoJet;
    curPseudoJet.reset_PtYPhiM(inParticles[i].pt,inParticles[i].eta,inParticles[i].phi,inParticles[i].m);
    curPseudoJet.set_user_index(inParticles[i].id);
    // fill vector of pseudojets
    _pfParticles.push_back(curPseudoJet);
    if(inParticles[i].id     > 3)  _genParticles.push_back( curPseudoJet);
    if(inParticles[i].id % 4 != 3) _pfchsParticles.push_back(curPseudoJet);
    // é
    if(inParticles[i].id % 4 == 2) _chargedPV.push_back(curPseudoJet);
    if(inParticles[i].id % 4 == 3) _chargedNoPV.push_back(curPseudoJet);
  }
}
puppiTMVAContainer::~puppiTMVAContainer() { }

double puppiTMVAContainer::goodVar(PseudoJet &iPart,std::vector<PseudoJet> &iParts, int iOpt) { 
  double Rsub = 0.3;
  double lPup = 0; 
  lPup = var_within_R(iOpt,iParts,iPart,Rsub);            
  if(iOpt == 6) lPup = lPup * iPart.pt()/pt_within_R(_pfParticles,iPart,Rsub);
  return lPup;
}
void puppiTMVAContainer::puppiForAll(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles) { 
  for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) { 
    double pVal = goodVar(iConstits[i0],iParticles,iOpt);
    _vals.push_back(pVal);
  }
}
double puppiTMVAContainer::compute(int iOpt,double iVal,double iMed,double iRMS,double iChi2Exp) { 
  if(iOpt == 1 && iVal < iMed) return 0;
  if(iOpt == 1 && iVal > iMed) return 1;
  double lVal = (iVal-iMed)/iRMS;
  int lId = 1; 
  lVal = lVal*fabs(lVal); 
  if(iChi2Exp != 0) lVal += iChi2Exp;
  if(iChi2Exp != 0) lId++;
  return  ROOT::Math::chisquared_cdf(lVal,lId);
}
void puppiTMVAContainer::computeMeanRMS() { 
  fMean = 0; 
  fRMS  = 0; 
  fMeanHEta = 0; 
  fRMSHEta  = 0; 
  int lNEvents = _recoParticles.size(); 
  std::vector<float> lValsPU;
  std::vector<float> lValsPUHEta;
  for(int i0 = 0; i0 < lNEvents; i0++) { 
    fVId      = float(_recoParticles[i0].vtxId);
    fPt       = _recoParticles[i0].pt;
    fEta      = _recoParticles[i0].eta;
    fTk       = _recoParticles[i0].trkChi2;
    fPVtx     = _recoParticles[i0].vtxChi2;
    fDepth    = _recoParticles[i0].depth;
    fTime     = _recoParticles[i0].time;
    fPFType   = float(_recoParticles[i0].pfType);
    fPuppi    = _vals[i0+0.*lNEvents];
    fPuppiO   = _vals[i0+1.*lNEvents];
    fPuppiLV  = _vals[i0+2.*lNEvents];
    fPuppiOLV = _vals[i0+3.*lNEvents];
    double pWeight = fReader     ->EvaluateMVA("BDT");
    pWeight += 1; pWeight *=0.5;
    _vals.push_back(pWeight);
    if(fabs(_recoParticles[i0].eta) > 2.3) lValsPUHEta.push_back(pWeight);
    if(fabs(_recoParticles[i0].eta) > 2.3) fMeanHEta += pWeight;
    if(fVId < 1) continue;
    fVId      = -1;
    pWeight = fReader     ->EvaluateMVA("BDT");
    pWeight += 1; pWeight *=0.5;
    fMean += pWeight;
    lValsPU.push_back(pWeight);
  }
  fMean/=lValsPU.size();
  for(unsigned int i0 = 0 ;i0 < lValsPU.size(); i0++) {fRMS += (lValsPU[i0]-fMean)*(lValsPU[i0]-fMean);}
  if(lValsPU.size() > 0)  fRMS/=lValsPU.size(); 
  fRMS = sqrt(fRMS);

  fMeanHEta/=lValsPUHEta.size();
  for(unsigned int i0 = 0 ;i0 < lValsPUHEta.size(); i0++) {fRMSHEta += (lValsPUHEta[i0]-fMeanHEta)*(lValsPUHEta[i0]-fMeanHEta);}
  if(lValsPUHEta.size() > 0)  fRMSHEta/=lValsPUHEta.size(); 
  fRMSHEta = sqrt(fRMSHEta);//*0.1;
}
std::vector<fastjet::PseudoJet> puppiTMVAContainer::puppiEvent     (int iOpt,double iQuant) { 
  std::vector<PseudoJet> particles;
  //Run through all compute mean and RMS
  _vals.resize(0);
  int lNEvents    = _pfParticles.size();
  puppiForAll(iOpt,_pfParticles,_pfParticles);
  puppiForAll(6   ,_pfParticles,_pfParticles);
  puppiForAll(iOpt,_pfParticles,_chargedPV  );
  puppiForAll(6   ,_pfParticles,_chargedPV  );
  computeMeanRMS();
  for(int i0 = 0; i0 < lNEvents; i0++) {
    fPt       = _recoParticles[i0].pt;
    fEta      = _recoParticles[i0].eta;
    fVId      = float(_recoParticles[i0].vtxId);
    double pChi2 = _recoParticles[i0].expChi2PU;
    if(_recoParticles[i0].pfType > 3) pChi2 = 0;
    double pWeight = _vals[i0+4.*lNEvents];
    if(fabs(_recoParticles[i0].eta)  < 2.6 ) pWeight =  compute(0,pWeight,fMean,fRMS,pChi2);
    if(fabs(_recoParticles[i0].eta)  > 2.6 ) pWeight =  compute(0,pWeight,fMeanHEta,fRMSHEta,pChi2);
    if(_recoParticles[i0].id % 4 == 3) pWeight = 0;
    if(_recoParticles[i0].id % 4 == 2) pWeight = 1;
    if(pWeight*fPt < 0.05) continue;
    PseudoJet curjet( pWeight*_pfParticles[i0].px(), pWeight*_pfParticles[i0].py(), pWeight*_pfParticles[i0].pz(), pWeight*_pfParticles[i0].e());
    curjet.set_user_index(_pfParticles[i0].user_index());
    particles.push_back(curjet);
  }
  return particles;
}
//FASTJET_END_NAMESPACE

