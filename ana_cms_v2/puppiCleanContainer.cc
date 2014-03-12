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
    fMin = 0.5;
    for (unsigned int i = 0; i < inParticles.size(); i++){
        fastjet::PseudoJet curPseudoJet;
        curPseudoJet.reset_PtYPhiM(inParticles[i].pt,inParticles[i].eta,inParticles[i].phi,inParticles[i].m);
        curPseudoJet.set_user_index(inParticles[i].id);
        // fill vector of pseudojets
        _pfParticles.push_back(curPseudoJet);

//        if(inParticles[i].id     > 3)  _genParticles.push_back( curPseudoJet);
//        if(inParticles[i].id % 4 != 3 && ((inParticles[i].pfType > 3 && curPseudoJet.pt() > fMin) || inParticles[i].pfType < 4)) _pfchsParticles.push_back(curPseudoJet);
//        // é
//        if(inParticles[i].id % 4 == 2) _chargedPV.push_back(curPseudoJet);
//        if(inParticles[i].id % 4 == 3) _chargedNoPV.push_back(curPseudoJet);

        if((inParticles[i].id == 0) && (inParticles[i].id == 2))  _genParticles.push_back( curPseudoJet);
        if(inParticles[i].id <= 2) _pfchsParticles.push_back(curPseudoJet);
        if(inParticles[i].id == 2) _chargedPV.push_back(curPseudoJet);
        if(inParticles[i].id == 3) _chargedNoPV.push_back(curPseudoJet);
        
    }
    
    puppiWeights_chLV.resize(0);
    puppiWeights_all.resize(0);
    alphas_chLV.resize(0);
    alphas_all.resize(0);
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
void puppiCleanContainer::getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,double iQuant,double iPtRMS,bool isForward) {
    std::vector<double> lValsPV;
    std::vector<double> lValsPU;
    for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) {
        double pVal = goodVar(iConstits[i0],iParticles,iOpt);
        _vals.push_back(pVal);
        if(iConstits[i0].pt() < iPtRMS) continue;//  && fabs(iConstits[i0].eta()) < 2.6*(iPtRMS == 0.5) ) continue;
        if( fabs(iConstits[i0].eta()) > 2.5  && isForward ) {lValsPU.push_back(pVal); continue;}
        if(isForward) continue;
        if(iConstits[i0].user_index() % 4 == 3) lValsPU.push_back(pVal);
        if(iConstits[i0].user_index() % 4 == 2) lValsPV.push_back(pVal);
    }
    std::sort (lValsPV.begin(),lValsPV.end());
    std::sort (lValsPU.begin(),lValsPU.end());
    double lMedPV = 0; if(lValsPV.size() > 0) lMedPV = lValsPV[int(lValsPV.size()*iQuant+0.5)];
    double lMedPU = 0; if(lValsPU.size() > 0) lMedPU = lValsPU[int(lValsPU.size()*iQuant+0.5)];
    double lRMSPU = 0;
    for(unsigned int i0 = 0 ;i0 < lValsPU.size(); i0++) {lRMSPU += (lValsPU[i0]-lMedPU)*(lValsPU[i0]-lMedPU);}
    double lRMSPV = 0;
    for(unsigned int i0 = 0 ;i0 < lValsPV.size(); i0++) {lRMSPV += (lValsPV[i0]-lMedPV)*(lValsPV[i0]-lMedPV);}
    if(lValsPV.size() > 0)  lRMSPV/=lValsPV.size();
    if(lValsPU.size() > 0)  lRMSPU/=lValsPU.size();
    fMed = lMedPU;
    fRMS = sqrt(lRMSPU);
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
std::vector<fastjet::PseudoJet> puppiCleanContainer::puppiEvent     (int iOpt,double iQuant) {
    std::vector<PseudoJet> particles;
    //Run through all compute mean and RMS
    _vals.resize(0);
    
    getRMSAvg(8,_pfParticles,_chargedPV,iQuant,0.5,false);
    double lMed0=fMed;
    double lRMS0=fRMS;
    int lNEvents    = _vals.size();
    
    getRMSAvg(8,_pfParticles,_pfParticles,iQuant,0.5,true);
    double lMed1=fMed;
    double lRMS1=fRMS;
    
    if(_isTuned) getRMSAvg(10,_pfParticles,_pfParticles,iQuant,0.5,true);
    double lMed2=fMed;
    double lRMS2=fRMS*0.5;
    for(int i0 = 0; i0 < lNEvents; i0++) {

        // fill alpha values
        alphas_chLV.push_back(_vals[i0]);
        alphas_all.push_back(_vals[i0+lNEvents]);
        
        puppiWeights_chLV.push_back( compute(0,_vals[i0]      ,lMed0,lRMS0,0)     );
        puppiWeights_all.push_back(  compute(0,_vals[i0+lNEvents]  ,lMed1,lRMS1,0) );
        
        double pWeight = 1;
        
        //?????
        double pChi2 = _recoParticles[i0].expChi2PU;
        if(!_isExperiment) pChi2 = 0;
        if(_recoParticles[i0].pfType > 3) pChi2 = 0;
        pChi2 = 0;
        
        if(fabs(_pfParticles[i0].eta()) < 2.5) pWeight *= compute(0,_vals[i0]            ,lMed0,lRMS0,pChi2);
        if(fabs(_pfParticles[i0].eta()) > 2.5) pWeight *= compute(0,_vals[i0+lNEvents]   ,lMed1,lRMS1,pChi2);
        //if(_isTuned) if(fabs(_pfParticles[i0].eta()) > 2.5) pWeight *= compute(0,_vals[i0+2.*lNEvents],lMed2,lRMS2,0);
        //CHS
        if(_pfParticles[i0].user_index() == 2 ) pWeight = 1;
        if(_pfParticles[i0].user_index() == 3 ) pWeight = 0;
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
