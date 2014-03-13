#include "puppiCleanContainer.hh"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "fastjet/internal/base.hh"
#include "TH2F.h"
#include "TMath.h"


using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


puppiCleanContainer::puppiCleanContainer(std::vector<RecoObj> inParticles,bool iExperiment,bool iTuned){
  _isExperiment = iExperiment;
  _isTuned      = iTuned;
  _recoParticles.resize(0); 
  _pfParticles.resize(0);
  _genParticles.resize(0);
  _pfchsParticles.resize(0);    
  _chargedPV.resize(0); 
  _chargedNoPV.resize(0);
  _vals.resize(0); 
  _recoParticles = inParticles;
  fMin = 1.0;
  for (unsigned int i = 0; i < inParticles.size(); i++){
    fastjet::PseudoJet curPseudoJet;
    curPseudoJet.reset_PtYPhiM(inParticles[i].pt,inParticles[i].eta,inParticles[i].phi,inParticles[i].m);
    curPseudoJet.set_user_index(inParticles[i].id);
    // fill vector of pseudojets
    _pfParticles.push_back(curPseudoJet);
    if(inParticles[i].id     > 3)  _genParticles.push_back( curPseudoJet);
    if(inParticles[i].id % 4 != 3 && ((inParticles[i].pfType > 3 && curPseudoJet.pt() > fMin) || inParticles[i].pfType < 4)) _pfchsParticles.push_back(curPseudoJet);
    // é
    if(inParticles[i].id % 4 == 2) _chargedPV.push_back(curPseudoJet);
    if(inParticles[i].id % 4 == 3) _chargedNoPV.push_back(curPseudoJet);
  }
}
puppiCleanContainer::~puppiCleanContainer(){}
double puppiCleanContainer::goodVar(PseudoJet &iPart,std::vector<PseudoJet> &iParts, int iOpt) { 
  double Rsub = 0.3; //For tracking using 0.2
  //double RLrg = 1.0;
  double lPup = 0; 
  if(iOpt == 10) return iPart.pt()+var_within_R(iOpt,iParts,iPart,Rsub);          
  lPup = var_within_R(iOpt,iParts,iPart,Rsub);            
  if(iOpt == 6) lPup = lPup * iPart.pt()/pt_within_R(_pfParticles,iPart,Rsub);
  return lPup;
}
void puppiCleanContainer::getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,double iQuant,double iPtRMS) { 
  std::vector<double> lValsPU;
  std::vector<double> lValsPUHEta;
  for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) { 
    double pVal = goodVar(iConstits[i0],iParticles,iOpt);
    _vals.push_back(pVal);
    if(iConstits[i0].pt() < iPtRMS) continue;//  && fabs(iConstits[i0].eta()) < 2.6*(iPtRMS == 0.5) ) continue;
    if( fabs(iConstits[i0].eta()) > 2.5   ) {lValsPUHEta.push_back(pVal); }
    if(iConstits[i0].user_index() % 4 == 3) lValsPU.push_back(pVal);
  }
  std::sort (lValsPU    .begin(),lValsPU    .end());   
  std::sort (lValsPUHEta.begin(),lValsPUHEta.end());  
  double lMedPU     = 0; if(lValsPU    .size() > 0) lMedPU     = lValsPU    [int(lValsPU    .size()*iQuant+0.5)];
  double lMedPUHEta = 0; if(lValsPUHEta.size() > 0) lMedPUHEta = lValsPUHEta[int(lValsPUHEta.size()*iQuant+0.5)]; 
  double lRMSPU     = 0; 
  for(unsigned int i0 = 0 ;i0 < lValsPU    .size(); i0++) {lRMSPU     += (lValsPU    [i0]-lMedPU)    *(lValsPU    [i0]-lMedPU);}
  double lRMSPUHEta = 0; 
  for(unsigned int i0 = 0 ;i0 < lValsPUHEta.size(); i0++) {lRMSPUHEta += (lValsPUHEta[i0]-lMedPUHEta)*(lValsPUHEta[i0]-lMedPUHEta);}
  if(lValsPU    .size() > 0)  lRMSPU    /=lValsPU.size(); 
  if(lValsPUHEta.size() > 0)  lRMSPUHEta/=lValsPUHEta.size(); 
  fMed     = lMedPU;  
  fRMS     = sqrt(lRMSPU);
  fMedHEta = lMedPUHEta;  
  fRMSHEta = sqrt(lRMSPUHEta);
}
double puppiCleanContainer::compute(int iOpt,double iVal,double iMed,double iRMS,double iChi2Exp) { 
  if(iOpt == 1 && iVal < iMed) return 0;
  if(iOpt == 1 && iVal > iMed) return 1;
  double lVal = (iVal-iMed)/iRMS;
  int lNDof = 1;
  if(iChi2Exp > 0) lVal += iChi2Exp;
  if(iChi2Exp > 0) lNDof++;
  return  ROOT::Math::chisquared_cdf(lVal*fabs(lVal),lNDof);
}
double puppiCleanContainer::compute2(int iOpt,double iVal,double iMed,double iRMS,double iChi2Exp,double iVal1,double iMed1,double iRMS1) { 
  if(iOpt == 1 && iVal < iMed) return 0;
  if(iOpt == 1 && iVal > iMed) return 1;
  double lVal  = (iVal -iMed) /iRMS;
  double lVal1 = (iVal1-iMed1)/iRMS1;
  int lNDof = 2;
  if(iChi2Exp > 0) lVal += iChi2Exp;
  if(iChi2Exp > 0) lNDof++;
  return  ROOT::Math::chisquared_cdf(lVal*fabs(lVal)+lVal1*fabs(lVal1),lNDof);
}
double puppiCleanContainer::compute3(int iOpt,double iVal,double iMed,double iRMS,double iChi2Exp,double iVal1,double iMed1,double iRMS1,double iVal2,double iMed2,double iRMS2) { 
  if(iOpt == 1 && iVal < iMed) return 0;
  if(iOpt == 1 && iVal > iMed) return 1;
  double lVal  = (iVal -iMed) /iRMS;
  double lVal1 = (iVal1-iMed1)/iRMS1;
  double lVal2 = (iVal2-iMed2)/iRMS2;
  int lNDof = 3;
  if(iChi2Exp > 0) lVal += iChi2Exp;
  if(iChi2Exp > 0) lNDof++;
  return  ROOT::Math::chisquared_cdf(lVal*fabs(lVal)+lVal1*fabs(lVal1)+lVal2*fabs(lVal2),lNDof);
}
double puppiCleanContainer::getChi2FromdZ(double iDZ) { 
  //We need to obtain prob of PU + (1-Prob of LV)
  // Prob(LV) = Gaus(dZ,sigma) where sigma = 1.5mm  (its really more like 1mm)
  double lProbLV = ROOT::Math::normal_cdf_c(fabs(iDZ),0.15)*2.; //*2 is to do it double sided
  double lProbPU = 1-lProbLV;
  if(lProbPU <= 0) lProbPU = 1e-16; 
  if(lProbPU >= 0) lProbPU = 1-1e-16; 
  double lChi2PU = TMath::ChisquareQuantile(lProbPU,1);
  return lChi2PU;
}
std::vector<fastjet::PseudoJet> puppiCleanContainer::puppiEvent     (int iOpt,double iQuant) { 
  std::vector<PseudoJet> particles;
  //Run through all compute mean and RMS
  _vals.resize(0);
  getRMSAvg(8,_pfParticles,_chargedPV,iQuant,0.5);
  double lMed0=fMed; 
  double lRMS0=fRMS;  
  int lNEvents    = _vals.size();
  
  getRMSAvg(8,_pfParticles,_pfParticles,iQuant,0.5);
  double lMed1=fMed; 
  double lRMS1=fRMS;
  double lMed1HEta=fMedHEta; 
  double lRMS1HEta=fRMSHEta;

  if(_isTuned) getRMSAvg(10,_pfParticles,_pfParticles,iQuant,0.5);
  double lMed2=fMed; 
  double lRMS2=fRMS;
  double lMed2HEta=fMedHEta; 
  double lRMS2HEta=fRMSHEta;
  for(int i0 = 0; i0 < lNEvents; i0++) {
    double pWeight = 1;
    double pChi2 = getChi2FromdZ(_recoParticles[i0].dZ);
    if(!_isExperiment) pChi2 = 0;
    if(_recoParticles[i0].pfType > 3) pChi2 = 0;
    if(!_isTuned && fabs(_pfParticles[i0].eta()) < 2.5) pWeight *= compute (0,_vals[i0],lMed0,lRMS0,pChi2);
    if( _isTuned && fabs(_pfParticles[i0].eta()) < 2.5) pWeight *= compute (0,_vals[i0],lMed0,lRMS0,pChi2);
    //if(_isTuned && fabs(_pfParticles[i0].eta()) < 2.5) pWeight *= compute (0,_vals[i0+lNEvents],lMed1,lRMS1,0);
    //if( _isTuned && fabs(_pfParticles[i0].eta()) < 2.5) pWeight *= compute2(0,_vals[i0]            ,lMed0,lRMS0,pChi2,_vals[i0+1.*lNEvents],lMed1,lRMS1);
    //if( _isTuned && fabs(_pfParticles[i0].eta()) < 2.5) pWeight *= compute3(0,_vals[i0]            ,lMed0,lRMS0,pChi2,_vals[i0+1.*lNEvents],lMed1,lRMS1,_vals[i0+2.*lNEvents],lMed2,lRMS2);
    
    if(!_isTuned) if(fabs(_pfParticles[i0].eta()) > 2.5) pWeight *= compute(0,_vals[i0+1.*lNEvents],lMed1HEta,lRMS1HEta,0);
    if(_isTuned)  if(fabs(_pfParticles[i0].eta()) > 2.5) pWeight *= compute(0,_vals[i0+1.*lNEvents],lMed1HEta,lRMS1HEta,0);//,lRMS1HEta,0,_vals[i0+2.*lNEvents],lMed2HEta,lRMS2HEta);
    if(_isTuned)  if(fabs(_pfParticles[i0].eta()) > 2.5) pWeight *= compute(0,_vals[i0+2.*lNEvents],lMed2HEta,lRMS2HEta,0);
    //CHS
    if(_pfParticles[i0].user_index() % 4 == 2 ) pWeight = 1;  
    if(_pfParticles[i0].user_index() % 4 == 3 ) pWeight = 0;
    //Basic Cuts
    if(pWeight*_pfParticles[i0].pt()  < 0.05) continue;
    if(pWeight*_pfParticles[i0].pt()  < fMin && _recoParticles[i0].pfType > 3 ) continue;
    //Produce
    PseudoJet curjet( pWeight*_pfParticles[i0].px(), pWeight*_pfParticles[i0].py(), pWeight*_pfParticles[i0].pz(), pWeight*_pfParticles[i0].e());
    curjet.set_user_index(_recoParticles[i0].id);
    particles.push_back(curjet);
  }
  return particles;
}
 //if(_pfParticles[i0].pt()  > 10  && _recoParticles[i0].pfType > 3 && fabs(_pfParticles[i0].eta()) > 2.5) pWeight = 1;
//FASTJET_END_NAMESPACE

