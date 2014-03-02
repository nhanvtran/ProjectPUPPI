#include "puppiContainer.hh"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "TH2F.h"
#include "TMath.h"
#include "fastjet/internal/base.hh"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


puppiContainer::puppiContainer(std::vector<fastjet::PseudoJet> inParticles, std::vector<int> isPU, std::vector<int> isCh,bool iDiscretize){
    
  _isPU = isPU;
  _isCh = isCh;
  
  _genParticles.resize(0);
  _pfchsParticles.resize(0);    
  _chargedPV.resize(0); 
  _chargedNoPV.resize(0);
  
  for (unsigned int i = 0; i < inParticles.size(); i++){
    // fill vector of pseudojets
    if (isPU[i] == 0){
      _genParticles.push_back( inParticles[i] );            
    }
    //if(inParticles[i].user_index() < 2 && iDiscretize) continue; //just charged
    _pfParticles.push_back(inParticles[i]);
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
  //if(iDiscretize) discretize(_pfParticles,inParticles);
}
void puppiContainer::discretize(std::vector<fastjet::PseudoJet> &discreteParticles,std::vector<fastjet::PseudoJet> &iParticles) {
  //Discretize in 0.2x0.2
  TH2F* h2d      = new TH2F( "h2d"     ,";#eta;#phi;e (GeV)",100, -5,5, 63,0,2*TMath::Pi() ); //100 63
  TH2F* h2dPU    = new TH2F( "h2dPU"   ,";#eta;#phi;e (GeV)",100, -5,5, 63,0,2*TMath::Pi() ); //100 63
  for (unsigned int i0 = 0; i0 < iParticles.size(); i0++){
    int curBin = h2d->FindBin(iParticles[i0].eta(),iParticles[i0].phi());
    if (iParticles[i0].user_index() == 0) h2d     ->SetBinContent( curBin, h2d     ->GetBinContent(curBin) + iParticles[i0].e() );        
    if (iParticles[i0].user_index() == 1) h2dPU   ->SetBinContent( curBin, h2dPU   ->GetBinContent(curBin) + iParticles[i0].e() );        
  }
  //For now take the bin center and not the weighted center
  for(int i0   = 0; i0 < h2d->GetNbinsX()+1; i0++) { 
    for(int i1 = 0; i1 < h2d->GetNbinsY()+1; i1++) { 
      double pEta   = h2d->GetXaxis()->GetBinCenter(i0);
      double pPhi   = h2d->GetYaxis()->GetBinCenter(i1);
      double pPt    = fabs(h2d  ->GetBinContent(i0,i1)*sin(atan(exp(-pEta))*2.));
      double pPtPU  = fabs(h2dPU->GetBinContent(i0,i1)*sin(atan(exp(-pEta))*2.));
      PseudoJet pJet;
      pJet.reset_PtYPhiM(pPt+pPtPU,pEta,pPhi);
      pJet.set_user_index(0);
      if(pPt < pPtPU ) pJet.set_user_index(1);
      if(pPt >  0 && pPtPU  > 0) pJet.set_user_index(pJet.user_index()+4); 
      if(pJet.pt() > 0.2) discreteParticles.push_back(pJet);
    }
  }
  delete h2d;
  delete h2dPU;
}

puppiContainer::~puppiContainer(){}

std::vector<fastjet::PseudoJet> puppiContainer::trimEvent(){
    
    //std::cout << "trimming event..." << std::endl;
    std::vector<PseudoJet> answer;
    answer.resize(0);
    // -- event trimming parameters -- 
    double Rjet = 1;
    double ptcut = 25;
    double Rsub = 0.4;
    double fcut = 0.3;	
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
double puppiContainer::goodVar(PseudoJet &iPart,std::vector<PseudoJet> &iParts, int iOpt) { 
  double Rsub = 0.3;
  double lPup = 0; 
  lPup = var_within_R(iOpt,iParts,iPart,Rsub);            
  if(iOpt > 5) lPup = lPup * iPart.pt()/pt_within_R(_pfParticles,iPart,Rsub);
  return lPup;
}


