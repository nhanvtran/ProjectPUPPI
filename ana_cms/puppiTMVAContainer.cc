#include "puppiTMVAContainer.hh"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "fastjet/internal/base.hh"
#include "TH2F.h"
#include "TMath.h"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

puppiTMVAContainer::puppiTMVAContainer(std::vector<fastjet::PseudoJet> &inParticles, std::vector<int> &isPU, std::vector<int> &isCh, std::vector<float> &iTk, std::vector<float> &iAVtx,std::vector<int> &iVId,std::string iWeight,bool iVtx,bool iDiscretize,bool iCHS) { 
  fCHS        = iCHS;
  fVtx        = iVtx;
  fDiscretize = iDiscretize;
  TMVA::Tools::Instance();
  fReader      = new TMVA::Reader( "!Color:!Silent" );    
  fPt       = 0; fReader->AddVariable("pt"                 , &fPt);
  fEta      = 0; fReader->AddVariable("eta"                , &fEta);
  fTk       = 0; fReader->AddVariable("tk"                 , &fTk);
  fPVtx     = 0; fReader->AddVariable("vtx"                , &fPVtx);
  fVId      = 0; fReader->AddVariable("vid"                , &fVId);
  fPuppi    = 0; fReader->AddVariable("ptodR"              , &fPuppi);
  fPuppiO   = 0; fReader->AddVariable("ptodRSO"            , &fPuppiO);
  fPuppiLV  = 0; if(fVtx) fReader->AddVariable("ptodR_lv"           , &fPuppiLV);
  fPuppiOLV = 0; if(fVtx) fReader->AddVariable("ptodRSO_lv"         , &fPuppiOLV);
  fReader     ->BookMVA("BDT",iWeight.c_str());
  refresh(inParticles,isPU,isCh,iTk,iAVtx,iVId); 
}
void puppiTMVAContainer::refresh(std::vector<fastjet::PseudoJet> &inParticles, std::vector<int> &isPU, std::vector<int> &isCh, std::vector<float> &iTk, std::vector<float> &iVtx,std::vector<int> &iVId) { 
  _pfParticles.resize(0);
  _genParticles.resize(0);
  _pfchsParticles.resize(0);    
  _chargedPV.resize(0); 
  _chargedNoPV.resize(0);
  _vals.resize(0); 
  _isPU = isPU;
  _isCh = isCh;
  _tk   = iTk;
  _vtx  = iVtx;
  _vid  = iVId;
  for (unsigned int i = 0; i < inParticles.size(); i++){
    // fill vector of pseudojets
    //if (fabs(inParticles[i].pt()) < 0.2) continue;        
    //if (fabs(inParticles[i].eta()) > 5)  continue;
    if (_isPU[i] == 0) _genParticles.push_back( inParticles[i] );            
    if (fabs(inParticles[i].user_index()) < 2 && fDiscretize)  continue;
    _pfParticles.push_back(inParticles[i]);
    if (!(_vid[i] > 0)) { // || (inParticles[i].pt() < 1 && fabs(inParticles[i].user_index()) < 2))) { 
      _pfchsParticles.push_back( inParticles[i] );
    }
    if (_vid[i] == 0) { 
      _chargedPV.push_back( inParticles[i] );
    }
    if (_vid[i] >  0) { 
      _chargedNoPV.push_back( inParticles[i] );
    }
  }
  //if(fDiscretize) discretize(_pfParticles   ,inParticles);
  //if(fDiscretize) discretize(_pfchsParticles,inParticles,true);
}
void puppiTMVAContainer::discretize(std::vector<fastjet::PseudoJet> &discreteParticles,std::vector<fastjet::PseudoJet> &iParticles,bool iPtCut) {
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
  return;
  std::sort (lValsPV.begin(),lValsPV.end());   
  std::sort (lValsPU.begin(),lValsPU.end());   
  double lMedPV = lValsPV[int(lValsPV.size()*iQuant+0.5)];
  double lMedPU = lValsPU[int(lValsPU.size()*iQuant+0.5)];
  std::cout << "NEvents 1 " << iConstits.size() << " --  " << _vals.size() << std::endl;
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
  int lNEvents    = _pfParticles.size();
  if(!fCHS) getRMSAvg(iOpt,_pfParticles,_pfParticles,_isPU,iQuant);
  if(!fCHS) getRMSAvg(6,_pfParticles,_pfParticles,_isPU,iQuant);
  if(!fCHS && fVtx) getRMSAvg(iOpt,_pfParticles,_chargedPV,_isPU,iQuant);
  if(!fCHS && fVtx) getRMSAvg(6   ,_pfParticles,_chargedPV,_isPU,iQuant);
  std::vector<double> lValsPU;
  float lMean = 0; 
  for(int i0 = 0; i0 < lNEvents; i0++) { 
    fVId      = float(_vid[i0]);
    if(fVId < 1) continue;
    fVId      = -1;
    fPt       = _pfParticles[i0].pt();
    fEta      = _pfParticles[i0].eta();
    fTk       = _tk[i0];
    fPVtx     = _vtx[i0];
    fPuppi    = _vals[i0+0.*lNEvents];
    fPuppiO   = _vals[i0+1.*lNEvents];
    fPuppiLV  = _vals[i0+2.*lNEvents];
    fPuppiOLV = _vals[i0+3.*lNEvents];
    double pWeight = fReader     ->EvaluateMVA("BDT");
    if(!fCHS) { pWeight += 1; pWeight *=0.5; }
    lMean += pWeight;
    lValsPU.push_back(pWeight);
  }
  lMean/=lValsPU.size();
  double lRMSPU = 0; 
  for(unsigned int i0 = 0 ;i0 < lValsPU.size(); i0++) {lRMSPU += (lValsPU[i0]-lMean)*(lValsPU[i0]-lMean);}
  if(lValsPU.size() > 0)  lRMSPU/=lValsPU.size(); 
  lRMSPU = sqrt(lRMSPU);
  cout << "Mean - " << lMean << " -- " << lRMSPU << endl;
  for(int i0 = 0; i0 < lNEvents; i0++) {
    fPt       = _pfParticles[i0].pt();
    fEta      = _pfParticles[i0].eta();
    fTk       = _tk[i0];
    fPVtx     = _vtx[i0];
    fVId      = float(_vid[i0]);
    if(fVId > 1) continue;
    double pWeight = 1;
    if(!fCHS) fPuppi    = _vals[i0+0.*lNEvents];
    if(!fCHS) fPuppiO   = _vals[i0+1.*lNEvents];
    if(!fCHS && fVtx) fPuppiLV  = _vals[i0+2.*lNEvents];
    if(!fCHS && fVtx) fPuppiOLV = _vals[i0+3.*lNEvents];
    if(!fCHS && fabs(_pfParticles[i0].eta()) < 12.5) pWeight = fReader     ->EvaluateMVA("BDT");
    if(!fCHS) { pWeight += 1; pWeight *=0.5; }
    if(fVId == 0) pWeight = 1;
    //if(fCHS) pWeight = 1;
    //if(_pfParticles[i0].user_index() == 2 ) pWeight = 1;  
    //if(_pfParticles[i0].user_index() == 3 ) pWeight = 0;
    //if(pWeight < 0.2) continue;
    if(fVId != 0) pWeight =  compute(0,pWeight,lMean,lRMSPU*2.);
    if(fVId > -1) pWeight *= fTk;
    if(fVId > -1) pWeight *= fPVtx;
    if(_pfParticles[i0].user_index()  < 2 && pWeight*fPt < 0.2) continue;
    PseudoJet curjet( pWeight*_pfParticles[i0].px(), pWeight*_pfParticles[i0].py(), pWeight*_pfParticles[i0].pz(), pWeight*_pfParticles[i0].e());
    curjet.set_user_index(_pfParticles[i0].user_index());
    particles.push_back(curjet);
  }
  return particles;
}
//FASTJET_END_NAMESPACE
