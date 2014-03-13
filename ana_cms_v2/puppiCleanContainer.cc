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
    fNeutralMinE  = 0.1;  //=> This can be tuned
    //Clear everything
    _recoParticles.resize(0);
    _pfParticles.resize(0);
    _genParticles.resize(0);
    _pfchsParticles.resize(0);
    _chargedPV.resize(0);
    _chargedNoPV.resize(0);
    _vals.resize(0);
    //Link to the RecoObjects
    _recoParticles = inParticles;
    for (unsigned int i = 0; i < inParticles.size(); i++){
        fastjet::PseudoJet curPseudoJet;
        curPseudoJet.reset_PtYPhiM(inParticles[i].pt,inParticles[i].eta,inParticles[i].phi,inParticles[i].m);
        curPseudoJet.set_user_index(inParticles[i].id);
        // fill vector of pseudojets for internal references
        _pfParticles.push_back(curPseudoJet);
        if((inParticles[i].id == 0) && (inParticles[i].id == 2))  _genParticles.push_back( curPseudoJet);
        if(inParticles[i].id <= 2) _pfchsParticles.push_back(curPseudoJet);                                 //Remove Charged particles associated to other vertex
        if(inParticles[i].id == 2) _chargedPV.push_back(curPseudoJet);                                      //Take Charged particles associated to PV
        if(inParticles[i].id == 3) _chargedNoPV.push_back(curPseudoJet);
        
    }
    puppiWeights_chLV.resize(0);
    puppiWeights_all.resize(0);
    alphas_chLV.resize(0);
    alphas_all.resize(0);
}
puppiCleanContainer::~puppiCleanContainer(){}

double puppiCleanContainer::goodVar(PseudoJet &iPart,std::vector<PseudoJet> &iParts, int iOpt) {
    double Rsub = 0.3;
    double lPup = 0;
    //Cryptic way of calling one of the Puppi Algos (there are indeed 10)
    if(iOpt == 10) return iPart.pt()+var_within_R(iOpt,iParts,iPart,Rsub);
    lPup = var_within_R(iOpt,iParts,iPart,Rsub);
    if(iOpt == 6) lPup = lPup * iPart.pt()/pt_within_R(_pfParticles,iPart,Rsub);
    return lPup;
}
//In fact takes the median no the average
void puppiCleanContainer::getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,double iQuant,double iPtRMS) { 
  std::vector<double> lValsPU;
  std::vector<double> lValsPUHEta;
  for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) { 
    double pVal = goodVar(iConstits[i0],iParticles,iOpt);
    _vals.push_back(pVal);
    if(iConstits[i0].pt() < iPtRMS) continue;
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
double puppiCleanContainer::getChi2FromdZ(double iDZ) { 
  //We need to obtain prob of PU + (1-Prob of LV)
  // Prob(LV) = Gaus(dZ,sigma) where sigma = 1.5mm  (its really more like 1mm)
  double lProbLV = ROOT::Math::normal_cdf_c(fabs(iDZ),0.2)*2.; //*2 is to do it double sided
  double lProbPU = 1-lProbLV;
  if(lProbPU <= 0) lProbPU = 1e-16;   //Quick Trick to through out infs
  if(lProbPU >= 0) lProbPU = 1-1e-16; //Ditto
  double lChi2PU = TMath::ChisquareQuantile(lProbPU,1);
  return lChi2PU;
}
std::vector<fastjet::PseudoJet> puppiCleanContainer::puppiEvent     (int iOpt,double iQuant) {
    std::vector<PseudoJet> particles;
    //Run through all compute mean and RMS
    _vals.resize(0);
    int lNEvents    = _recoParticles.size();
    
    //Log2 Puppi with Leading Vertex Particles Only
    getRMSAvg(8,_pfParticles,_chargedPV,iQuant,0.5);
    double lMed0=fMed;
    double lRMS0=fRMS;
    
    //Log2 Puppi with all Particles
    getRMSAvg(8,_pfParticles,_pfParticles,iQuant,0.5);
    double lMed1=fMed;
    double lRMS1=fRMS;
    double lMed1HEta=fMedHEta; 
    double lRMS1HEta=fRMSHEta;

    //Sum pT in a cone 
    if(_isTuned) getRMSAvg(10,_pfParticles,_pfParticles,iQuant,0.5);
    double lMed2=fMed;
    double lRMS2=fRMS;
    double lMed2HEta=fMedHEta; 
    double lRMS2HEta=fRMSHEta;
    for(int i0 = 0; i0 < lNEvents; i0++) {
        // fill alpha values
        alphas_chLV.push_back(_vals[i0]);
        alphas_all.push_back(_vals[i0+lNEvents]);
        // fill the p-values
        puppiWeights_chLV.push_back( compute(0,_vals[i0]      ,lMed0,lRMS0,0)     );
        puppiWeights_all.push_back(  compute(0,_vals[i0+lNEvents]  ,lMed1,lRMS1,0) );
        //Now compute the puppi Weight
        double pWeight = 1;
	double pChi2   = 0;
        if(_isExperiment) {
	  //Compute an Experimental Puppi Weight with delta Z info (very simple example)
	  pChi2 = getChi2FromdZ(_recoParticles[i0].dZ);
	  //Now make sure Neutrals are not set
	  if(_recoParticles[i0].pfType > 3) pChi2 = 0;
	}
	//Basic Puppi
        if(fabs(_pfParticles[i0].eta()) < 2.5) pWeight *= compute(0,_vals[i0]            ,lMed0,lRMS0,pChi2);
        if(fabs(_pfParticles[i0].eta()) > 2.5) pWeight *= compute(0,_vals[i0+lNEvents]   ,lMed1HEta,lRMS1HEta,pChi2);
        if(_isTuned) if(fabs(_pfParticles[i0].eta()) > 2.5) pWeight *= compute(0,_vals[i0+2.*lNEvents],lMed2HEta,lRMS2HEta,0);
        //CHS
        if(_pfParticles[i0].user_index() == 2 ) pWeight = 1;
        if(_pfParticles[i0].user_index() == 3 ) pWeight = 0;
        //Basic Cuts
        if(pWeight*_pfParticles[i0].pt()  < 0.05) continue;  //==> Elminate the low pT stuff 
        if(pWeight*_pfParticles[i0].E()   < fNeutralMinE && _recoParticles[i0].pfType > 3 ) continue;  //threshold cut on the neutral E
        
        //Produce
        PseudoJet curjet( pWeight*_pfParticles[i0].px(), pWeight*_pfParticles[i0].py(), pWeight*_pfParticles[i0].pz(), pWeight*_pfParticles[i0].e());
        curjet.set_user_index(_recoParticles[i0].id);
        particles.push_back(curjet);

    }
    return particles;
}


