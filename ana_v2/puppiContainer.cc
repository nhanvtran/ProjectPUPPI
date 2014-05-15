#include "puppiContainer.hh"

#include "fastjet/internal/base.hh"
#include "Math/ProbFunc.h"
#include "TH2F.h"
#include "TMath.h"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


puppiContainer::puppiContainer(std::vector<fastjet::PseudoJet> inParticles, std::vector<int> isPU, std::vector<int> isCh, bool iDiscretize){
    
    //_pfParticles = inParticles;
    _isPU = isPU;
    _isCh = isCh;
    _isPFCHS.resize(0);
    
    _pfParticles.resize(0);
    _genParticles.resize(0);
    _pfchsParticles.resize(0); 
    _pfchsParticlesToUse.resize(0);    
       
    _neutrals.resize(0);
    _chargedLV.resize(0);    
    _chargedPU.resize(0);    
    
    std::vector< fastjet::PseudoJet > neutralParts;
    neutralParts.clear();
    
    double etaTracker = 2.5;
    double userIndex_isPuCh = -1;
    
    int ctr_00 = 0;
    int ctr_01 = 0;
    int ctr_10 = 0;
    int ctr_11 = 0;            
    
    for (unsigned int i = 0; i < inParticles.size(); i++){
        float weightGEN = 0.;
        float weightPFCHS = 0.;
    
        if (isPU[i] == 0 && isCh[i] == 0) ctr_00++;
        if (isPU[i] == 0 && isCh[i] == 1) ctr_01++;
        if (isPU[i] == 1 && isCh[i] == 0) ctr_10++;
        if (isPU[i] == 1 && isCh[i] == 1) ctr_11++;
                                        
        // fill vector of pseudojets
        if (fabs(inParticles[i].eta()) < 5 && inParticles[i].pt() > 0.0){        
            
            PseudoJet curjet(inParticles[i].px(), inParticles[i].py(), inParticles[i].pz(), inParticles[i].e());
            curjet.set_user_index( -99 );
            
            if (isPU[i] == 0){ // is GEN
                PseudoJet curjetGEN( inParticles[i].px(), inParticles[i].py(), inParticles[i].pz(), inParticles[i].e());
                curjetGEN.set_user_index( 0 );
                _genParticles.push_back( curjetGEN );
            }
            
            if ((isPU[i] == 0) || (isPU[i] == 1 && isCh[i] == 0 && fabs(inParticles[i].eta()) < etaTracker) || (isPU[i] == 1 && fabs(inParticles[i].eta()) > etaTracker)){
                weightPFCHS = 1.;
            }

            if (weightPFCHS > 0){
                if (isCh[i] == 0 || fabs(inParticles[i].eta()) > etaTracker){
                    if (isPU[i] == 0) userIndex_isPuCh = 0;
                    else userIndex_isPuCh = 1;
                }
                else{
                    userIndex_isPuCh = 2;                    
                }
            }
            else{
                userIndex_isPuCh = 3;      
            }
            
            //std::cout << "isPU[i] = " << isPU[i] << ", isCh[i] = " << isCh[i] << ", userIndex_isPuCh = " << userIndex_isPuCh << ", inParticles[i].eta()) = " << inParticles[i].eta() << std::endl;
            curjet.set_user_index( userIndex_isPuCh );
            
            if (userIndex_isPuCh == 0 || userIndex_isPuCh == 1) neutralParts.push_back( curjet );
            else if (userIndex_isPuCh == 2) _chargedLV.push_back( curjet );                    
            else if (userIndex_isPuCh == 3)_chargedPU.push_back( curjet );
            else { std::cout << "extra particles?" << std::endl; }
            
            _isPFCHS.push_back( weightPFCHS );
            
        }
    }
    
    puppiWeights_pfchs.resize(0);
    puppiWeights_chLV.resize(0);
    puppiWeights_all.resize(0);
    alphas_chLV.resize(0);
    alphas_all.resize(0);
    cleansedWeights.resize(0);    
    
    // discretize neutrals? 
    bool bDiscretize = iDiscretize;
    std::vector< fastjet::PseudoJet > neutralParts_discrete;
        
    if (bDiscretize){
        discretize(neutralParts_discrete,neutralParts);
        _neutrals = neutralParts_discrete;
    }
    else _neutrals = neutralParts;

//    std::cout << "neutralParts_discrete = " << neutralParts_discrete.size() << std::endl;
//    std::cout << "neutralParts = " << neutralParts.size() << std::endl;    
//    double totalEnergy = 0;
//    double totalEnergy_original = 0;
//    for (unsigned int i = 0; i < neutralParts.size(); i++){
//        totalEnergy_original += neutralParts[i].pt();
//    }
//    for (unsigned int i = 0; i < neutralParts_discrete.size(); i++){
//        totalEnergy += neutralParts_discrete[i].pt();
//    }
//    std::cout << "totalEnergy = " << totalEnergy << ", totalEnergy_original = " << totalEnergy_original << std::endl;
    
    // ===========================
    // set pf particles
//    std::cout << " _neutrals.size() " << _neutrals.size() << std::endl;           
//    std::cout << " _chargedLV.size() " << _chargedLV.size() << std::endl;           
//    std::cout << " _chargedPU.size() " << _chargedPU.size() << std::endl;   
    double ptcutOnNeutrals = 0.0; 
    for (unsigned int i = 0; i < _neutrals.size(); i++){
        _pfParticles.push_back(_neutrals[i]);
        _pfchsParticles.push_back(_neutrals[i]);  
        if (_neutrals[i].pt() > ptcutOnNeutrals){
            _pfchsParticlesToUse.push_back(_neutrals[i]);  
        }
    }
    for (unsigned int i = 0; i < _chargedLV.size(); i++){
        _pfParticles.push_back(_chargedLV[i]);
        _pfchsParticles.push_back(_chargedLV[i]); 
        _pfchsParticlesToUse.push_back(_chargedLV[i]);                
                       
    }
    for (unsigned int i = 0; i < _chargedPU.size(); i++){
        _pfParticles.push_back(_chargedPU[i]);        
    }    
    // ===========================    
    
//    _pfParticles.clear();
//    _pfParticles = inParticles;

}

puppiContainer::~puppiContainer(){}

void puppiContainer::discretize(std::vector<fastjet::PseudoJet> &discreteParticles,std::vector<fastjet::PseudoJet> &iParticles) {
  bool lGrid = false;
  //Discretize in 0.1x0.1
//  TH2F* h2d      = new TH2F( "h2d"     ,";#eta;#phi;e (GeV)",2000, -5,5, 1250,0,2*TMath::Pi() ); //100 63
//  TH2F* h2dPU    = new TH2F( "h2dPU"   ,";#eta;#phi;e (GeV)",2000, -5,5, 1250,0,2*TMath::Pi() ); //100 63
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
    int curBin = h2d->FindBin(iParticles[i0].eta(),iParticles[i0].phi());
    if (fabs(iParticles[i0].pt())  < 0.01) continue;        
    if (fabs(iParticles[i0].eta()) > 5)   continue;        
    if (lGrid  && iParticles[i0].user_index() ==  0) h2d     ->SetBinContent( curBin, h2d     ->GetBinContent(curBin) + iParticles[i0].e() );        
    if (lGrid  && iParticles[i0].user_index() ==  1) h2dPU   ->SetBinContent( curBin, h2dPU   ->GetBinContent(curBin) + iParticles[i0].e() );        
    if (!lGrid && iParticles[i0].user_index() ==  0) lJets  [curBin] += iParticles[i0];
    if (!lGrid && iParticles[i0].user_index() ==  1) lJetsPU[curBin] += iParticles[i0];
  }
  //For now sumit up
  if(!lGrid) { 
    for(unsigned int i0   = 0; i0 < lJets.size(); i0++) { 
      PseudoJet pJet;
      if(lJets  [i0].pt() > 0) pJet   += lJets  [i0]; 
      if(lJetsPU[i0].pt() > 0) pJet   += lJetsPU[i0];
      pJet.set_user_index(0);
      if(lJets[i0].pt() < lJetsPU[i0].pt() ) pJet.set_user_index(1);
      if(pJet.pt() > 0.1) discreteParticles.push_back(pJet);
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
	if(pJet.pt() > 0.1) discreteParticles.push_back(pJet);
      }
    }
  }
  delete h2d;
  delete h2dPU;
}

std::vector<fastjet::PseudoJet> puppiContainer::trimEvent(double vRjet, double vptcut, double vRsub, double vfcut ){
    
    //std::cout << "trimming event..." << std::endl;
    std::vector<PseudoJet> answer;
    answer.resize(0);
    
    // -- event trimming parameters -- 
    double Rjet = vRjet;
    double ptcut = vptcut;
    double Rsub = vRsub;
    double fcut = vfcut;	
    
    for(unsigned int i=0; i<_pfParticles.size(); i++){
        double weight0 = 0.0;        
        
        // --- event trimming ---
        double pt_Rjet = pt_within_R(_pfParticles,_pfParticles[i],Rjet);
        double pt_Rsub = 0;
        if (pt_Rjet >= ptcut){
            
            pt_Rsub = pt_within_R(_pfParticles,_pfParticles[i],Rsub);            
            if (pt_Rsub / pt_Rjet > fcut) weight0 = 1.0;
        }
        if (weight0 > 0){
            PseudoJet curjet( weight0*_pfParticles[i].px(), weight0*_pfParticles[i].py(), weight0*_pfParticles[i].pz(), weight0*_pfParticles[i].e());
            answer.push_back( curjet );
        }
    }
    
    return answer;
}

std::vector<fastjet::PseudoJet> puppiContainer::cleanseEvent( double Rsub ){
    
    std::vector<PseudoJet> answer;
    answer.resize(0);
    cleansedWeights.resize(0);    
    
    // -- event cleansing parameters -- 
    double Rsub_cleansing = Rsub;
    
    std::vector<PseudoJet> charged_PU;
    std::vector<PseudoJet> charged_LV;
    for(unsigned int i = 0; i<_pfParticles.size(); i++){
        if(_isPU[i] == 1 && _isCh[i] == 1) charged_PU.push_back(_pfParticles[i]);
        if(_isPU[i] == 0 && _isCh[i] == 1) charged_LV.push_back(_pfParticles[i]);
    }    
    
    for(unsigned int i=0; i<_pfParticles.size(); i++){
        double weight1 = 0.0;
        
        // --- event cleansing ---        
        double pt_R_CPU = pt_within_R(charged_PU,_pfParticles[i],Rsub_cleansing);
        double pt_R_CLV = pt_within_R(charged_LV,_pfParticles[i],Rsub_cleansing);        
        
        if (pt_R_CLV + pt_R_CPU > 0) weight1 = pt_R_CLV / (pt_R_CLV + pt_R_CPU );
        if (weight1 > 0){
            PseudoJet curjet( weight1*_pfParticles[i].px(), weight1*_pfParticles[i].py(), weight1*_pfParticles[i].pz(), weight1*_pfParticles[i].e());
            answer.push_back( curjet );
        }
        cleansedWeights.push_back( weight1 );    
    }
    
    return answer;
}

std::vector<fastjet::PseudoJet> puppiContainer::puppiEvent_V1( double Rsub, double exponent ){
    
    std::vector<PseudoJet> answer;
    answer.resize(0);
    puppiWeights_all.resize(0);
    
    // -- puppi cleansing parameters -- 
    double Rsub_cleansing = Rsub;
    
    std::vector<PseudoJet> charged_PU;
    std::vector<PseudoJet> charged_LV;
    for(unsigned int i = 0; i<_pfParticles.size(); i++){
        if(_isPU[i] == 1 && _isCh[i] == 1) charged_PU.push_back(_pfParticles[i]);
        if(_isPU[i] == 0 && _isCh[i] == 1) charged_LV.push_back(_pfParticles[i]);
    }    
    
    for(unsigned int i=0; i<_pfParticles.size(); i++){
        double weight1 = 0.0;
        
        // --- puppi cleansing ---        
        double pt_R_CPU = ktWeight_within_R(charged_PU,_pfParticles[i],Rsub_cleansing,exponent);
        double pt_R_CLV = ktWeight_within_R(charged_LV,_pfParticles[i],Rsub_cleansing,exponent);        
        
        if (pt_R_CLV + pt_R_CPU > 0) weight1 = pt_R_CLV / (pt_R_CLV + pt_R_CPU );
        if (weight1 > 0){
            PseudoJet curjet( weight1*_pfParticles[i].px(), weight1*_pfParticles[i].py(), weight1*_pfParticles[i].pz(), weight1*_pfParticles[i].e());
            answer.push_back( curjet );
        }
        puppiWeights_all.push_back( weight1 );    
    }
    
    return answer;
}

// =================================================
// =================================================
// =================================================
// =================================================
// =================================================

std::vector<fastjet::PseudoJet> puppiContainer::puppiEvent(int iPU, double iQuant){
    std::vector<PseudoJet> particles;
    //Run through all compute mean and RMS
    _vals.resize(0);
    
    double R0 = 0.3;
    double R1 = 0.3;
    
    // the chi2 2dof version
    getRMSAvg(13,_pfParticles,_chargedLV,_isPU,iQuant,0.5,R0);
    //getRMSAvg(4,_pfParticles,_chargedLV,_isPU,iQuant,0.5,R0);
    double lMed0=fMed;
    double lRMS0=fRMS;
    int lNEvents  = _vals.size();
    
    getRMSAvg(13,_pfParticles,_pfParticles,_isPU,0.5,0.2,R1);
    //double lMed1=fMedHEta;
    //double lRMS1=fRMSHEta;
    double lMed1=fMed;
    double lRMS1=fRMS;
    
    float wptCutC = 0.5;
    float wptCutF = 1.0;
//    if ((iPU < 30)&&(iPU > 15)){ wptCutC = 0.15; wptCutF = 0.25; }
//    if ((iPU < 45)&&(iPU > 30)){ wptCutC = 0.15; wptCutF = 0.4; }
//    if ((iPU < 70)&&(iPU >= 45)){ wptCutC = 0.3; wptCutF = 0.65; }
//    if ((iPU < 100)&&(iPU >= 70)){ wptCutC = 0.4; wptCutF = 0.8; }
//    if ((iPU < 120)&&(iPU >= 100)){ wptCutC = 0.65; wptCutF = 1.25;}
//    if ((iPU < 150)&&(iPU >= 120)){ wptCutC = 0.8; wptCutF = 1.5; }
    
    //a functional form, hard-coded for now
    wptCutC = 0.66667e-2*( (float) iPU ) + 0.1;
    wptCutF = 1.05e-2*( (float) iPU ) + 0.2;
    
    for(int i0 = 0; i0 < lNEvents; i0++) {
        double pWeight = 1;
        
        // fill alpha values
        alphas_chLV.push_back(_vals[i0]);
        alphas_all.push_back(_vals[i0+lNEvents]);
        
        puppiWeights_chLV.push_back( compute(0,_vals[i0]           ,lMed0,lRMS0) );
        puppiWeights_all.push_back(  compute(0,_vals[i0+lNEvents]  ,lMed1,lRMS1) );
        
        if(fabs(_pfParticles[i0].eta()) < 2.5){
            //pWeight *= compute2dof(_vals[i0],lMed0,lRMS0,_vals[i0+lNEvents],lMed1,lRMS1);
            pWeight *= compute(0,_vals[i0],lMed0,lRMS0);
        }
        if(fabs(_pfParticles[i0].eta()) > 2.5){
            pWeight *= compute(0,_vals[i0+lNEvents]   ,lMed1,lRMS1);
        }
        //if(fabs(_pfParticles[i0].eta()) > 2.5) pWeight *= compute(0,_vals[i0+lNEvents*2.],lMed2,lRMS2);
        if(_pfParticles[i0].user_index() == 2 ) pWeight = 1;
        if(_pfParticles[i0].user_index() == 3 ) pWeight = 0;
        if(_pfParticles[i0].user_index()  < 2 && pWeight*_pfParticles[i0].pt() < wptCutC && fabs(_pfParticles[i0].eta()) < 2.5) continue;
        if(_pfParticles[i0].user_index()  < 2 && pWeight*_pfParticles[i0].pt() < wptCutF && fabs(_pfParticles[i0].eta()) > 2.5) continue;
        if(pWeight < 0.1) continue;
        
        PseudoJet curjet( pWeight*_pfParticles[i0].px(), pWeight*_pfParticles[i0].py(), pWeight*_pfParticles[i0].pz(), pWeight*_pfParticles[i0].e());
        curjet.set_user_index(_pfParticles[i0].user_index());
        particles.push_back(curjet);
    }
    return particles;
}
double puppiContainer::goodVar(PseudoJet &iPart,std::vector<PseudoJet> &iParts, int iOpt, double R0) {
    double Rsub = R0; //For tracking using 0.2
                       //double RLrg = 1.0;
    double lPup = 0; 
    lPup = var_within_R(iOpt,iParts,iPart,Rsub);            
    if(iOpt == 6) lPup = lPup * iPart.pt()/pt_within_R(_pfParticles,iPart,Rsub);
    return lPup;
}
void puppiContainer::getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,std::vector<int> &iIsPU,double iQuant,double iPtRMS, double R0, bool countZeros) {

    ////std::vector<double> lValsPV;
    std::vector<double> lValsPU;
    std::vector<double> lValsPUHEta;
    
    std::vector<double> curVals;
    for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) {
        double pVal = goodVar(iConstits[i0],iParticles,iOpt,R0);
        curVals.push_back(pVal);
        if(iConstits[i0].pt() < iPtRMS) continue;
        if(pVal == 0) continue;
        if(iConstits[i0].user_index()  == 3) lValsPU.push_back(pVal);
        ////if(iConstits[i0].user_index()  == 2) lValsPV.push_back(pVal);
        if( fabs(iConstits[i0].eta()) > 2.5   ) {lValsPUHEta.push_back(pVal); }
    }

    std::sort (lValsPU.begin(),lValsPU.end());
    std::sort (lValsPUHEta.begin(),lValsPUHEta.end());

    double lMedPU = 0; if(lValsPU.size() > 0) lMedPU = lValsPU[int(lValsPU.size()*iQuant+0.5)];
    double lMedPUHEta = 0; if(lValsPUHEta.size() > 0) lMedPUHEta = lValsPUHEta[int(lValsPUHEta.size()*iQuant+0.5)];
    
    // now add back the 0 particles to the vectorr to compute the RMS
    for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) {
//        //std::cout << "curVals[i0] = " << curVals[i0] << std::endl;
//        if (curVals[i0] != 0) _vals.push_back(curVals[i0]);
//        else _vals.push_back(lMedPU);
        _vals.push_back(curVals[i0]);
    }
    
    // now compute RMS
    double lRMSPU = 0;
    int lRMSPUctr = 0;
    for(unsigned int i0 = 0 ;i0 < lValsPU.size(); i0++) {
        if ((lValsPU[i0]-lMedPU) > 0) continue;
        lRMSPU += (lValsPU[i0]-lMedPU)*(lValsPU[i0]-lMedPU);
        lRMSPUctr++;
    }
    
    double lRMSPUHEta = 0;
    int lRMSPUHEtactr = 0;
    for(unsigned int i0 = 0 ;i0 < lValsPUHEta.size(); i0++) {
        //if ((lValsPUHEta[i0]-lMedPUHEta) > 0) continue;
        lRMSPUHEta += (lValsPUHEta[i0]-lMedPUHEta)*(lValsPUHEta[i0]-lMedPUHEta);
        lRMSPUHEtactr++;
    }
    
    ////if(lValsPV.size() > 0)  lRMSPV/=lValsPV.size();
    if(lValsPU.size() > 0)  lRMSPU/=lRMSPUctr;
    if(lValsPUHEta.size() > 0)  lRMSPUHEta/=lRMSPUHEtactr;
    
    fMed = lMedPU;
    fRMS = sqrt(lRMSPU);
    fMedHEta = lMedPUHEta;
    fRMSHEta = sqrt(lRMSPUHEta);
    
}
double puppiContainer::compute(int iOpt,double iVal,double iMed,double iRMS) { 
    if(iOpt == 1 && iVal < iMed) return 0;
    if(iOpt == 1 && iVal > iMed) return 1;
    double lVal = (iVal-iMed)/iRMS;
    //return lVal*fabs(lVal);
    //return TMath::Erf(lVal);//*fabs(lVal);
    return  ROOT::Math::chisquared_cdf(lVal*fabs(lVal),1.);
}

double puppiContainer::compute2dof(double iVal,double iMed,double iRMS,double iVal2,double iMed2,double iRMS2) {
    double lVal = (iVal-iMed)/iRMS;
    double lVal2 = (iVal2-iMed2)/iRMS2;
    return  ROOT::Math::chisquared_cdf(lVal*fabs(lVal) + lVal2*fabs(lVal2),2.);
}
// =================================================

double puppiContainer::pt_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R){
    
    fastjet::Selector sel = fastjet::SelectorCircle(R);
    sel.set_reference(centre);
    vector<PseudoJet> near_particles = sel(particles);
    double answer = 0.0;
    //std::cout << "near particles (pt) = " << near_particles.size() << std::endl;
    
    for(unsigned int i=0; i<near_particles.size(); i++){
        answer += near_particles[i].pt();
    }
    return(answer);
}

double puppiContainer::ktWeight_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R, double exponent){
    
    fastjet::Selector sel = fastjet::SelectorCircle(R);
    sel.set_reference(centre);
    vector<PseudoJet> near_particles = sel(particles);
    
    double sumWeight = 0.0;
    //std::cout << "near particles (kt) = " << near_particles.size() << std::endl;
    for(unsigned int i=0; i<near_particles.size(); i++){
        double deta = (centre.eta() - near_particles[i].eta());
        double dphi = (centre.phi() - near_particles[i].phi());
        double drij2 = deta*deta + dphi*dphi;
        double drij = sqrt(drij2);
        if (drij2 > 0){
            //float ktFact = min( pow(near_particles[i].pt(),2), pow(centre.pt(),2) )*drij2;
            float weight = min( pow(near_particles[i].pt(),exponent),pow(centre.pt(),exponent) ) / pow(drij,exponent);
            //std::cout << "ktFact = " << ktFact << std::endl;
            sumWeight += weight;
        }
    }
    return(sumWeight);
}

double puppiContainer::getAverage(const vector<double> & particles){
    //std::cout << "get mean" << std::endl;
    double ltotal = 0.;
    for(int i0 = 0 ;i0 < particles.size(); i0++) {ltotal += particles[i0];}
    return ltotal/particles.size();
}

double puppiContainer::getRMS(const vector<double> & particles){
    
    double average = getAverage( particles );
    double lRMS = 0;
    for(int i0 = 0 ;i0 < particles.size(); i0++) {lRMS += (particles[i0]-average)*(particles[i0]-average);}
    return lRMS/particles.size(); 
}


//FASTJET_END_NAMESPACE
