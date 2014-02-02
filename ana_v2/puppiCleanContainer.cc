#include "puppiCleanContainer.hh"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "fastjet/internal/base.hh"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


puppiCleanContainer::puppiCleanContainer(std::vector<fastjet::PseudoJet> inParticles, std::vector<int> isPU, std::vector<int> isCh){
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

puppiCleanContainer::~puppiCleanContainer(){}

std::vector<fastjet::PseudoJet> puppiCleanContainer::trimEvent(){
    
    //std::cout << "trimming event..." << std::endl;
    std::vector<PseudoJet> answer;
    answer.resize(0);
    
    // -- event trimming parameters -- 
    double Rjet = 1;
    double ptcut = 25;
    double Rsub = 0.4;
    double fcut = 0.3;	
    fMed = 0; 
    fRMS = 0; 
    for(unsigned int i=0; i<_pfParticles.size(); i++){
        double weight0 = 0.0;        
        
        // --- event trimming ---
        double pt_Rjet = pt_within_R(_pfParticles,_pfParticles[i],Rjet);
        double pt_Rsub = 0;
        if (pt_Rjet >= ptcut){
            
            pt_Rsub = pt_within_R(_pfParticles,_pfParticles[i],Rsub);            
            if (pt_Rsub / pt_Rjet > fcut) weight0 = 1.0;
        }
        if (weight0 != 0.0){
	  PseudoJet curjet( weight0*_pfParticles[i].px(), weight0*_pfParticles[i].py(), weight0*_pfParticles[i].pz(), weight0*_pfParticles[i].e());
	  answer.push_back( curjet );
        }
    }
    
    return answer;
}
double puppiCleanContainer::goodVar(PseudoJet &iPart,std::vector<PseudoJet> &iParts, int iOpt) { 
  double Rsub = 0.3;
  double lPup = 0; 
  lPup = var_within_R(iOpt,iParts,iPart,Rsub);            
  if(iOpt == 6) lPup = lPup * iPart.pt()/pt_within_R(_pfParticles,iPart,Rsub);
  return lPup;
}
void puppiCleanContainer::getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,std::vector<int> &iIsPU,double iQuant) { 
  std::vector<double> lValsPV;
  std::vector<double> lValsPU;
  for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) { 
    double pVal = goodVar(iConstits[i0],iParticles,iOpt);
    if( iIsPU[iConstits[i0].user_index()]) lValsPU.push_back(pVal);
    if(!iIsPU[iConstits[i0].user_index()]) lValsPV.push_back(pVal);
    _vals.push_back(pVal);
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
double puppiCleanContainer::compute(int iOpt,double iVal,double iMed,double iRMS) { 
  if(iOpt == 1 && iVal < iMed) return 0;
  if(iOpt == 1 && iVal > iMed) return 1;
  double lVal = (iVal-iMed)/iRMS;
  return  ROOT::Math::chisquared_cdf(lVal*lVal,1.);
}
std::vector<fastjet::PseudoJet> puppiCleanContainer::puppiEvent     (int iOpt,double iQuant) { 
  std::vector<PseudoJet> particles;
  //Run through all compute mean and RMS
  _vals.resize(0);
  getRMSAvg(iOpt,_pfParticles,_pfParticles,_isPU,iQuant);
  double lMed0=fMed; 
  double lRMS0=fRMS;  

  int lNEvents    = _vals.size();
  if(iOpt == 7) getRMSAvg(6,_pfParticles,_pfParticles,_isPU,iQuant);
  double lMed1=fMed; 
  double lRMS1=fRMS;  

  if(iOpt == 7) getRMSAvg(iOpt,_pfParticles,_chargedPV,_isPU,iQuant);
  double lMed2=fMed; 
  double lRMS2=fRMS;  

  for(int i0 = 0; i0 < lNEvents; i0++) {
    double pWeight = compute(1,_vals[i0]                    ,lMed0,lRMS0);
    if(iOpt == 7) pWeight *= compute(1,_vals[i0+lNEvents]   ,lMed1,lRMS1);
    if(iOpt == 7) pWeight *= compute(1,_vals[i0+2.*lNEvents],lMed2,lRMS2);
    if(pWeight < 0.1) continue;
    PseudoJet curjet( pWeight*_pfParticles[i0].px(), pWeight*_pfParticles[i0].py(), pWeight*_pfParticles[i0].pz(), pWeight*_pfParticles[i0].e());
    particles.push_back(curjet);
  }
  return particles;
}
//FASTJET_END_NAMESPACE

