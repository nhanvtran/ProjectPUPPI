#include "puppiTMVAContainer.hh"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "fastjet/internal/base.hh"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


puppiTMVAContainer::puppiTMVAContainer(std::vector<fastjet::PseudoJet> inParticles, std::vector<int> isPU, std::vector<int> isCh,std::string iWeight,bool iVtx) { 
  fVtx = iVtx;
  TMVA::Tools::Instance();
  fReader = new TMVA::Reader( "!Color:!Silent" );    
  fPt       = 0; fReader->AddVariable("pt"                 , &fPt);
  fPuppi    = 0; fReader->AddVariable("ptodR"              , &fPuppi);
  fPuppiO   = 0; fReader->AddVariable("ptodRSO"            , &fPuppiO);
  fPuppiLV  = 0; if(fVtx) fReader->AddVariable("ptodR_lv"           , &fPuppiLV);
  fPuppiOLV = 0; if(fVtx) fReader->AddVariable("ptodRSO_lv"         , &fPuppiOLV);
  fReader->BookMVA("BDT",iWeight.c_str());
  _isPU = isPU;
  _isCh = isCh;
  _pfParticles.resize(0);
  _genParticles.resize(0);
  _pfchsParticles.resize(0);    
  _chargedPV.resize(0); 
  _chargedNoPV.resize(0);
  _vals.resize(0); 
  
  for (unsigned int i = 0; i < inParticles.size(); i++){
    // fill vector of pseudojets
    if (fabs(inParticles[i].pt()) < 0.2) continue;        
    if (fabs(inParticles[i].eta()) < 5){        
      _pfParticles.push_back(inParticles[i]);
     if (isPU[i] == 0){
	_genParticles.push_back( inParticles[i] );            
      }
      if ((isPU[i] == 0) || (isPU[i] == 1 && isCh[i] == 0 && fabs(inParticles[i].eta()) < 2.5) || (isPU[i] == 1 && fabs(inParticles[i].eta()) > 2.5)){
	_pfchsParticles.push_back( inParticles[i] );
      }
      if ((isPU[i] == 0) && isCh[i] == 1 && fabs(inParticles[i].eta()) < 12.5){
	_chargedPV.push_back( inParticles[i] );
      }
      if ((isPU[i] == 1) && isCh[i] == 1  && fabs(inParticles[i].eta()) < 12.5){
	_chargedNoPV.push_back( inParticles[i] );
      }
    }
  }
}
void puppiTMVAContainer::refresh(std::vector<fastjet::PseudoJet> inParticles, std::vector<int> isPU, std::vector<int> isCh) { 
  _isPU = isPU;
  _isCh = isCh;
  _pfParticles.resize(0);
  _genParticles.resize(0);
  _pfchsParticles.resize(0);    
  _chargedPV.resize(0); 
  _chargedNoPV.resize(0);
  _vals.resize(0); 
  
  for (unsigned int i = 0; i < inParticles.size(); i++){
    // fill vector of pseudojets
    if (fabs(inParticles[i].pt()) < 0.2) continue;        
    if (fabs(inParticles[i].eta()) < 5){        
      _pfParticles.push_back(inParticles[i]);
     if (isPU[i] == 0){
	_genParticles.push_back( inParticles[i] );            
      }
      if ((isPU[i] == 0) || (isPU[i] == 1 && isCh[i] == 0 && fabs(inParticles[i].eta()) < 2.5) || (isPU[i] == 1 && fabs(inParticles[i].eta()) > 2.5)){
	_pfchsParticles.push_back( inParticles[i] );
      }
      if ((isPU[i] == 0) && isCh[i] == 1 && fabs(inParticles[i].eta()) < 12.5){
	_chargedPV.push_back( inParticles[i] );
      }
      if ((isPU[i] == 1) && isCh[i] == 1  && fabs(inParticles[i].eta()) < 12.5){
	_chargedNoPV.push_back( inParticles[i] );
      }
    }
  }
}
puppiTMVAContainer::~puppiTMVAContainer(){}

double puppiTMVAContainer::goodVar(PseudoJet &iPart,std::vector<PseudoJet> &iParts, int iOpt) { 
  double Rsub = 0.3;
  double lPup = 0; 
  lPup = var_within_R(iOpt,iParts,iPart,Rsub);            
  if(iOpt == 6) lPup = lPup * iPart.pt()/pt_within_R(_pfParticles,iPart,Rsub);
  return lPup;
}
void puppiTMVAContainer::getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,std::vector<int> &iIsPU,double iQuant) { 
  std::vector<double> lValsPV;
  std::vector<double> lValsPU;
  for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) { 
    double pVal = goodVar(iConstits[i0],iParticles,iOpt);
    _vals.push_back(pVal);
    if(iConstits[i0].pt() < 0.5) continue;
    if(iConstits[i0].user_index() % 2 == 1) lValsPU.push_back(pVal);
    if(iConstits[i0].user_index() % 2 == 0) lValsPV.push_back(pVal); 
    //if( iIsPU[iConstits[i0].user_index()]) lValsPU.push_back(pVal);
    //if(!iIsPU[iConstits[i0].user_index()]) lValsPV.push_back(pVal);
  }
  std::sort (lValsPV.begin(),lValsPV.end());   
  std::sort (lValsPU.begin(),lValsPU.end());   
  double lMedPV = lValsPV[int(lValsPV.size()*iQuant+0.5)];
  double lMedPU = lValsPU[int(lValsPU.size()*iQuant+0.5)];
  double lRMSPU = 0; 
  for(unsigned int i0 = 0 ;i0 < lValsPU.size(); i0++) {lRMSPU += (lValsPU[i0]-lMedPU)*(lValsPU[i0]-lMedPU);}
  double lRMSPV = 0; 
  for(unsigned int i0 = 0 ;i0 < lValsPV.size(); i0++) {lRMSPV += (lValsPV[i0]-lMedPV)*(lValsPV[i0]-lMedPV);}
  if(lValsPV.size() > 0)  lRMSPV/=lValsPV.size(); 
  if(lValsPU.size() > 0)  lRMSPU/=lValsPU.size(); 
  //iVals[0] = lMedPV;  
  //iVals[1] = sqrt(lRMSPV);
  fMed = lMedPU;  
  fRMS = sqrt(lRMSPU);
}
double puppiTMVAContainer::compute(int iOpt,double iVal,double iMed,double iRMS) { 
  if(iOpt == 1 && iVal < iMed) return 0;
  if(iOpt == 1 && iVal > iMed) return 1;
  double lVal = (iVal-iMed)/iRMS;
  return  ROOT::Math::chisquared_cdf(lVal*lVal,1.);
}
std::vector<fastjet::PseudoJet> puppiTMVAContainer::puppiEvent     (int iOpt,double iQuant) { 
  std::vector<PseudoJet> particles;
  //Run through all compute mean and RMS
  _vals.resize(0);
  getRMSAvg(iOpt,_pfParticles,_pfParticles,_isPU,iQuant);
  int lNEvents    = _vals.size();
  getRMSAvg(6,_pfParticles,_pfParticles,_isPU,iQuant);
  if(fVtx) getRMSAvg(iOpt,_pfParticles,_chargedPV,_isPU,iQuant);
  if(fVtx) getRMSAvg(6   ,_pfParticles,_chargedPV,_isPU,iQuant);

  for(int i0 = 0; i0 < lNEvents; i0++) {
    fPt       = _pfParticles[i0].pt();
    fPuppi    = _vals[i0+0.*lNEvents];
    fPuppiO   = _vals[i0+1.*lNEvents];
    if(fVtx) fPuppiLV  = _vals[i0+2.*lNEvents];
    if(fVtx) fPuppiOLV = _vals[i0+3.*lNEvents];
    double pWeight = fReader->EvaluateMVA("BDT");
    pWeight += 1; pWeight *=0.5;
    if(_pfParticles[i0].user_index() == 2) pWeight = 1;  
    if(_pfParticles[i0].user_index() == 3) pWeight = 0;
    if(_pfParticles[i0].user_index() <  2 && fPt < 1.5) pWeight = 0; 
    //if(pWeight < 0.5) pWeight = 0; 
    if(pWeight < 0.1) continue;
    pWeight = 1;
    PseudoJet curjet( pWeight*_pfParticles[i0].px(), pWeight*_pfParticles[i0].py(), pWeight*_pfParticles[i0].pz(), pWeight*_pfParticles[i0].e());
    particles.push_back(curjet);
  }
  return particles;
}
//FASTJET_END_NAMESPACE

