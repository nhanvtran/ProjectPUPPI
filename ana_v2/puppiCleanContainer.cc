#include "puppiCleanContainer.hh"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "fastjet/internal/base.hh"
#include "TH2F.h"
#include "TMath.h"


using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


puppiCleanContainer::puppiCleanContainer(std::vector<fastjet::PseudoJet> inParticles, std::vector<int> isPU, std::vector<int> isCh,bool iDiscretize){
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
    if (fabs(inParticles[i].eta()) > 5)  continue;
    if (_isPU[i] == 0) _genParticles.push_back( inParticles[i] );            
    if (fabs(inParticles[i].user_index()) < 2 && iDiscretize)  continue;
    _pfParticles.push_back(inParticles[i]);
    if (inParticles[i].user_index() != 3) { 
      _pfchsParticles.push_back( inParticles[i] );
    }
    if (inParticles[i].user_index() == 2) { 
      _chargedPV.push_back( inParticles[i] );
    }
    if (inParticles[i].user_index() == 3) { 
      _chargedNoPV.push_back( inParticles[i] );
    }
  }
  if(iDiscretize) discretize(_pfParticles   ,inParticles);
  if(iDiscretize) discretize(_pfchsParticles,inParticles,true);
}
void puppiCleanContainer::discretize(std::vector<fastjet::PseudoJet> &discreteParticles,std::vector<fastjet::PseudoJet> &iParticles,bool iPtCut) {
  bool lGrid = false;
  //Discretize in 0.1x0.1
  TH2F* h2d      = new TH2F( "h2d"     ,";#eta;#phi;e (GeV)",100, -5,5, 63,0,2*TMath::Pi() ); //100 63
  TH2F* h2dPU    = new TH2F( "h2dPU"   ,";#eta;#phi;e (GeV)",100, -5,5, 63,0,2*TMath::Pi() ); //100 63
  std::vector<PseudoJet> lJets; 
  std::vector<PseudoJet> lJetsPU; 
  for(int i0   = 0; i0 < h2d->GetNbinsX()+1; i0++) { 
    for(int i1 = 0; i1 < h2d->GetNbinsY()+1; i1++) { 
      PseudoJet pJet(0,0,0,0);
      PseudoJet pJetPU(0,0,0,0);
      lJets.push_back(pJet);
      lJetsPU.push_back(pJetPU);
    }
  }
  for (unsigned int i0 = 0; i0 < iParticles.size(); i0++){
    if (fabs(iParticles[i0].pt())  < 0.2) continue;        
    if (fabs(iParticles[i0].eta()) > 5  ) continue;        
    int curBin = h2d->FindBin(iParticles[i0].eta(),iParticles[i0].phi());
    if (lGrid  && iParticles[i0].user_index() == 0) h2d     ->SetBinContent( curBin, h2d     ->GetBinContent(curBin) + iParticles[i0].e() );        
    if (lGrid  && iParticles[i0].user_index() == 1) h2dPU   ->SetBinContent( curBin, h2dPU   ->GetBinContent(curBin) + iParticles[i0].e() );        
    if (!lGrid && iParticles[i0].user_index() == 0) lJets  [curBin] += iParticles[i0];
    if (!lGrid && iParticles[i0].user_index() == 1) lJetsPU[curBin] += iParticles[i0];
  }
  //For now sumit up
  if(!lGrid) { 
    for(unsigned int i0   = 0; i0 < lJets.size(); i0++) { 
      PseudoJet pJet;
      if(lJets  [i0].pt() > 0) pJet   += lJets  [i0]; 
      if(lJetsPU[i0].pt() > 0) pJet   += lJetsPU[i0];
      pJet.set_user_index(0);
      if(lJets[i0].pt() < lJetsPU[i0].pt() ) pJet.set_user_index(1);
      if(pJet.pt() < 0.2) continue; 
      if(pJet.pt() < 1.0 && iPtCut) continue; 
      discreteParticles.push_back(pJet);
    }
  }
  if(lGrid) { 
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
	if(pJet.pt() > 0.2) discreteParticles.push_back(pJet);
      }
    }
  }
  delete h2d;
  delete h2dPU;
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
  double Rsub = 0.3; //For tracking using 0.2
  //double RLrg = 1.0;
  double lPup = 0; 
  lPup = var_within_R(iOpt,iParts,iPart,Rsub);            
  if(iOpt == 6) lPup = lPup * iPart.pt()/pt_within_R(_pfParticles,iPart,Rsub);
  return lPup;
}
void puppiCleanContainer::getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,std::vector<int> &iIsPU,double iQuant,double iPtRMS) { 
  std::vector<double> lValsPV;
  std::vector<double> lValsPU;
  for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) { 
    double pVal = goodVar(iConstits[i0],iParticles,iOpt);
    _vals.push_back(pVal);
    if(iConstits[i0].pt() < iPtRMS) continue;
    if(iConstits[i0].user_index()  == 3) lValsPU.push_back(pVal);
    if(iConstits[i0].user_index()  == 2) lValsPV.push_back(pVal);
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
  //iVals[0] = lMedPV;  
  //iVals[1] = sqrt(lRMSPV);
  fMed = lMedPU;  
  fRMS = sqrt(lRMSPU);
}
double puppiCleanContainer::compute(int iOpt,double iVal,double iMed,double iRMS) { 
  if(iOpt == 1 && iVal < iMed) return 0;
  if(iOpt == 1 && iVal > iMed) return 1;
  double lVal = (iVal-iMed)/iRMS;
  //return lVal*fabs(lVal);
  //return TMath::Erf(lVal);//*fabs(lVal);
  return  ROOT::Math::chisquared_cdf(lVal*fabs(lVal),1.);
}
std::vector<fastjet::PseudoJet> puppiCleanContainer::puppiEvent     (int iOpt,double iQuant) { 
  std::vector<PseudoJet> particles;
  //Run through all compute mean and RMS
  _vals.resize(0);
  
  getRMSAvg(7,_pfParticles,_chargedPV,_isPU,iQuant,0.5);
  //getRMSAvg(6,_pfParticles,_pfParticles,_isPU,0.5);
  double lMed0=fMed; 
  double lRMS0=fRMS;  
  int lNEvents    = _vals.size();
  //if(iOpt == 7) getRMSAvg(6,_pfParticles,_pfParticles,_isPU,0.1,0.5);
  
  getRMSAvg(7,_pfParticles,_pfParticles,_isPU,0.5,0.2);
  double lMed1=fMed; 
  double lRMS1=fRMS;  

  //if(iOpt == 7) getRMSAvg(6,_pfParticles,_pfParticles,_isPU,iQuant);
  //double lMed2=fMed; 
  //double lRMS2=fRMS;  

  for(int i0 = 0; i0 < lNEvents; i0++) {
    double pWeight = 1;
    if(fabs(_pfParticles[i0].eta()) < 2.5) pWeight *= compute(0,_vals[i0]            ,lMed0,lRMS0);
    if(fabs(_pfParticles[i0].eta()) > 2.5) pWeight *= compute(0,_vals[i0+lNEvents]   ,lMed1,lRMS1);
    //if(iOpt == 7) pWeight *= compute(0,_vals[i0+lNEvents*2.],lMed2,lRMS2);
    //pWeight = ROOT::Math::chisquared_cdf(pWeight,2.);
    if(_pfParticles[i0].user_index() == 2 ) pWeight = 1;  
    if(_pfParticles[i0].user_index() == 3 ) pWeight = 0;
    if(_pfParticles[i0].user_index()  < 2 && pWeight*_pfParticles[i0].pt() < 1.0) continue;
    if(pWeight < 0.1) continue;
    PseudoJet curjet( pWeight*_pfParticles[i0].px(), pWeight*_pfParticles[i0].py(), pWeight*_pfParticles[i0].pz(), pWeight*_pfParticles[i0].e());
    particles.push_back(curjet);
  }
  return particles;
}
//FASTJET_END_NAMESPACE

