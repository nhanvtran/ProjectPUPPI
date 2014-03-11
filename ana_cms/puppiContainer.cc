#include "puppiContainer.hh"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "TH2F.h"
#include "TMath.h"
#include "fastjet/internal/base.hh"
#include "fastjet/PseudoJet.hh"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


puppiContainer::puppiContainer(std::vector<RecoObj> inParticles){
  _genParticles  .resize(0);
  _pfParticles   .resize(0);    
  _pfchsParticles.resize(0);    
  _chargedPV     .resize(0); 
  _chargedNoPV   .resize(0);
  for (unsigned int i = 0; i < inParticles.size(); i++){
    fastjet::PseudoJet curPseudoJet;
    curPseudoJet.reset_PtYPhiM(inParticles[i].pt,inParticles[i].eta,inParticles[i].phi,inParticles[i].m);
     // fill vector of pseudojets
    _pfParticles.push_back(curPseudoJet);
    if(inParticles[i].id > 4)  _genParticles.push_back( curPseudoJet);
    if(inParticles[i].id % 4 != 3) _pfchsParticles.push_back(curPseudoJet);
    if(inParticles[i].id % 4 == 2) _chargedPV.push_back(curPseudoJet);
    if(inParticles[i].id % 4 == 3) _chargedNoPV.push_back(curPseudoJet);
  }
}
puppiContainer::~puppiContainer(){}
double puppiContainer::goodVar(int iId,std::vector<PseudoJet> &iParts, int iOpt) { 
  double Rsub = 0.3;
  double lPup = 0; 
  lPup = var_within_R(iOpt,iParts,_pfParticles[iId],Rsub);            
  if(iOpt > 5) lPup = lPup * _pfParticles[iId].pt()/pt_within_R(_pfParticles,_pfParticles[iId],Rsub);
  return lPup;
}

