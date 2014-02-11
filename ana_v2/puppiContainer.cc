#include "puppiContainer.hh"

#include "fastjet/internal/base.hh"
#include "Math/ProbFunc.h"
#include "TH2F.h"
#include "TMath.h"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


puppiContainer::puppiContainer(std::vector<fastjet::PseudoJet> inParticles, std::vector<int> isPU, std::vector<int> isCh){
    
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
    cleansedWeights.resize(0);    
    
    // discretize neutrals? 
    bool bDiscretize = true;
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

// ================================================= Nhan's old code

//std::vector<fastjet::PseudoJet> puppiContainer::puppiEvent(int iPU, double iQuant){
//    
//    double Rsub = 0.3;///Using 0.3 radius
//    std::vector<PseudoJet> answer;
//    answer.resize(0);
//
//    std::vector<double> lVals;
//    std::vector<double> lValsPV;
//    std::vector<double> lValsAll;
//    std::vector<double> lValsAPV;
//    std::vector<double> lValsAPU;
//    std::vector<double> lValsAAll;    
//    double lAvgPV = 0; 
//    double lAvgPU = 0; 
//    double lAvgAll = 0;
//    double lRMSPU = 0; 
//    double lRMSPV = 0;      
//    double lRMSAll = 0; 
//        
//    //adding the pt_i/Sum(pt_j) metric
//    std::vector<double> lVals_wpi;
//    std::vector<double> lValsPV_wpi;
//    std::vector<double> lValsAll_wpi;
//    std::vector<double> lValsAPV_wpi;
//    std::vector<double> lValsAPU_wpi;
//    std::vector<double> lValsAAll_wpi;    
//    double lAvgPV_wpi = 0; 
//    double lAvgPU_wpi = 0; 
//    double lAvgAll_wpi = 0;         
//    double lRMSPU_wpi = 0; 
//    double lRMSPV_wpi = 0;      
//    double lRMSAll_wpi = 0; 
//        
//    for(unsigned int i=0; i<_pfParticles.size(); i++){
//        
//        double puppi_Rsub = puppi_within_R(_pfchsParticles,_pfParticles[i],Rsub,false);      
//        double puppi_RsubPV = puppi_within_R(_chargedLV,_pfParticles[i],Rsub,false);   
//        double puppi_RsubAll = puppi_within_R(_pfParticles,_pfParticles[i],Rsub,false);   
//        
//        // downweighting for soft particles near hard particles
//        double pt_Rsub = pt_within_R(_pfParticles,_pfParticles[i],Rsub);
//        double weightI = _pfParticles[i].pt()/pt_Rsub;
////        
////        puppi_Rsub*=weightI;
////        puppi_RsubPV*=weightI;        
//        
//        lVals.push_back(puppi_Rsub);
//        lValsPV.push_back(puppi_RsubPV);
//        lValsAll.push_back(puppi_RsubAll);
//        
//        lVals_wpi.push_back(puppi_Rsub*weightI);
//        lValsPV_wpi.push_back(puppi_RsubPV*weightI);
//        lValsAll_wpi.push_back(puppi_RsubAll*weightI);        
//                
//        //Was applying a pt cut
//        if(_pfParticles[i].pt() > 0.2 && _pfParticles[i].user_index() == 3) { 
//            lValsAPV.push_back(puppi_RsubPV);
//            lValsAPU.push_back(puppi_Rsub);
//            lValsAAll.push_back(puppi_RsubAll);
//            
//            lValsAPV_wpi.push_back(puppi_RsubPV*weightI);
//            lValsAPU_wpi.push_back(puppi_Rsub*weightI);
//            lValsAAll_wpi.push_back(puppi_RsubAll*weightI);
//        }
//        
//    }
//    
//    //    //shift to make positive-definite
//    //    for(int i0 = 0; i0 < lVals.size(); i0++) { 
//    //        lVals[i0] -= min_lVals;
//    //        lValsPV[i0] -= min_lValsPV;
//    //    }
//    
//    //Median
//    std::sort (lValsAPV.begin(),lValsAPV.end());  
//    std::sort (lValsAPU.begin(),lValsAPU.end());  
//    std::sort (lValsAAll.begin(),lValsAAll.end());  
//    lAvgPV = lValsAPV[int(lValsAPV.size()/2.+0.5)];
//    lAvgPU = lValsAPU[int(lValsAPU.size()/2.+0.5)];
//    lAvgAll = lValsAAll[int(lValsAAll.size()/2.+0.5)];
//    //RMS
//    for(int i0 = 0 ;i0 < lValsAPU.size(); i0++) {lRMSPU += (lValsAPU[i0]-lAvgPU)*(lValsAPU[i0]-lAvgPU);}
//    for(int i0 = 0 ;i0 < lValsAPV.size(); i0++) {lRMSPV += (lValsAPV[i0]-lAvgPV)*(lValsAPV[i0]-lAvgPV);}
//    for(int i0 = 0 ;i0 < lValsAAll.size(); i0++) {lRMSAll += (lValsAAll[i0]-lAvgAll)*(lValsAAll[i0]-lAvgAll);}
//    if(lValsAPV.size() > 0) lRMSPV/=lValsAPV.size(); 
//    if(lValsAPU.size() > 0) lRMSPU/=lValsAPU.size(); 
//    if(lValsAAll.size() > 0) lRMSAll/=lValsAAll.size(); 
//
//    //Median
//    std::sort (lValsAPV_wpi.begin(),lValsAPV_wpi.end());  
//    std::sort (lValsAPU_wpi.begin(),lValsAPU_wpi.end());  
//    std::sort (lValsAAll_wpi.begin(),lValsAAll_wpi.end());  
//    lAvgPV_wpi = lValsAPV_wpi[int(lValsAPV_wpi.size()/2.+0.5)];
//    lAvgPU_wpi = lValsAPU_wpi[int(lValsAPU_wpi.size()/2.+0.5)];
//    lAvgAll_wpi = lValsAAll_wpi[int(lValsAAll_wpi.size()/2.+0.5)];
//    //RMS
//    for(int i0 = 0 ;i0 < lValsAPU_wpi.size(); i0++) {lRMSPU_wpi += (lValsAPU_wpi[i0]-lAvgPU_wpi)*(lValsAPU_wpi[i0]-lAvgPU_wpi);}
//    for(int i0 = 0 ;i0 < lValsAPV_wpi.size(); i0++) {lRMSPV_wpi += (lValsAPV_wpi[i0]-lAvgPV_wpi)*(lValsAPV_wpi[i0]-lAvgPV_wpi);}
//    for(int i0 = 0 ;i0 < lValsAAll_wpi.size(); i0++) {lRMSAll_wpi += (lValsAAll_wpi[i0]-lAvgAll_wpi)*(lValsAAll_wpi[i0]-lAvgAll_wpi);}
//    if(lValsAPV_wpi.size() > 0) lRMSPV_wpi/=lValsAPV_wpi.size(); 
//    if(lValsAPU_wpi.size() > 0) lRMSPU_wpi/=lValsAPU_wpi.size(); 
//    if(lValsAAll_wpi.size() > 0) lRMSAll_wpi/=lValsAAll_wpi.size(); 
//    
//    //    std::cout << "min_lVals = " << min_lVals << std::endl;
//    //    std::cout << "min_lValsPV = " << min_lValsPV << std::endl;
//    
//    //Particles
//    for(int i0 = 0; i0 < lVals.size(); i0++) { 
//        
//        double pWeight = lVals[i0];
//        double pWeight1 = lValsPV[i0];
//        double pWeight2 = lValsAll[i0];
//        
//        double chiWeight = (lVals[i0]-lAvgPU)/sqrt(lRMSPU)*(fabs(lVals[i0]-lAvgPU)/sqrt(lRMSPU)); 
//        double chiWeight1 = (lValsPV[i0]-lAvgPV)/sqrt(lRMSPV)*(fabs(lValsPV[i0]-lAvgPV)/sqrt(lRMSPV)); 
//        double chiWeight2 = (lValsAll[i0]-lAvgAll)/sqrt(lRMSAll)*(fabs(lValsAll[i0]-lAvgAll)/sqrt(lRMSAll)); 
//        pWeight = ROOT::Math::chisquared_cdf(chiWeight,1.);
//        pWeight1 = ROOT::Math::chisquared_cdf(chiWeight1,1.);
//        pWeight2 = ROOT::Math::chisquared_cdf(chiWeight2,1.);
//
//        double chiWeight_wpi = (lVals_wpi[i0]-lAvgPU_wpi)/sqrt(lRMSPU_wpi)*(fabs(lVals_wpi[i0]-lAvgPU_wpi)/sqrt(lRMSPU_wpi)); 
//        double chiWeight1_wpi = (lValsPV_wpi[i0]-lAvgPV_wpi)/sqrt(lRMSPV_wpi)*(fabs(lValsPV_wpi[i0]-lAvgPV_wpi)/sqrt(lRMSPV_wpi)); 
//        double chiWeight2_wpi = (lValsAll_wpi[i0]-lAvgAll_wpi)/sqrt(lRMSAll_wpi)*(fabs(lValsAll_wpi[i0]-lAvgAll_wpi)/sqrt(lRMSAll_wpi)); 
//        double pWeight_wpi = ROOT::Math::chisquared_cdf(chiWeight_wpi,1.);
//        double pWeight1_wpi = ROOT::Math::chisquared_cdf(chiWeight1_wpi,1.);
//        double pWeight2_wpi = ROOT::Math::chisquared_cdf(chiWeight2_wpi,1.);
//        
////        pWeight *= pWeight_wpi;
////        pWeight1 *= pWeight1_wpi;
////        pWeight2 *= pWeight2_wpi;
//                
//        //std::cout << "pWeight2 = " << pWeight2 << std::endl;
//        
//        puppiWeights_pfchs.push_back(pWeight);
//        puppiWeights_chLV.push_back(pWeight1);    
//        puppiWeights_all.push_back(pWeight2);    
//        
//        double weightToUse = 1.0;
//        if ( fabs(_pfParticles[i0].eta()) > 2.5 ) weightToUse *= pWeight;
//        if ( fabs(_pfParticles[i0].eta()) < 2.5 ) weightToUse *= pWeight1;
//        if(_pfParticles[i0].user_index() == 2 ) weightToUse = 1; 
//        if(_pfParticles[i0].user_index() == 3 ) weightToUse = 0;
//        if(_pfParticles[i0].user_index() < 2 && weightToUse*_pfParticles[i0].pt() < 1.0) continue;
//        if(weightToUse < 0.1) continue;
//        PseudoJet curjet(weightToUse*_pfParticles[i0].px(), 
//                         weightToUse*_pfParticles[i0].py(), 
//                         weightToUse*_pfParticles[i0].pz(), 
//                         weightToUse*_pfParticles[i0].e());
//        curjet.set_user_index(_pfParticles[i0].user_index());            
//        answer.push_back( curjet );      
//    }
//    
//    return answer;
//}

// =================================================

std::vector<fastjet::PseudoJet> puppiContainer::puppiEvent(int iPU, double iQuant){
    std::vector<PseudoJet> particles;
    //Run through all compute mean and RMS
    _vals.resize(0);
    
    getRMSAvg(7,_pfParticles,_chargedLV,_isPU,iQuant,0.5);
    //getRMSAvg(6,_pfParticles,_pfParticles,_isPU,0.5);
    double lMed0=fMed; 
    double lRMS0=fRMS; 
    int lNEvents  = _vals.size();
    //if(iOpt == 7) getRMSAvg(6,_pfParticles,_pfParticles,_isPU,0.1,0.5);
    
    getRMSAvg(7,_pfParticles,_pfParticles,_isPU,0.5,0.2);
    double lMed1=fMed; 
    double lRMS1=fRMS; 
    
    //if(iOpt == 7) getRMSAvg(6,_pfParticles,_pfParticles,_isPU,iQuant);
    //double lMed2=fMed; 
    //double lRMS2=fRMS; 
    
    for(int i0 = 0; i0 < lNEvents; i0++) {
        double pWeight = 1;
        if(fabs(_pfParticles[i0].eta()) < 2.5) pWeight *= compute(0,_vals[i0]      ,lMed0,lRMS0);
        if(fabs(_pfParticles[i0].eta()) > 2.5) pWeight *= compute(0,_vals[i0+lNEvents]  ,lMed1,lRMS1);
        //if(iOpt == 7) pWeight *= compute(0,_vals[i0+lNEvents*2.],lMed2,lRMS2);
        //pWeight = ROOT::Math::chisquared_cdf(pWeight,2.);
        if(_pfParticles[i0].user_index() == 2 ) pWeight = 1; 
        if(_pfParticles[i0].user_index() == 3 ) pWeight = 0;
        if(_pfParticles[i0].user_index() < 2 && pWeight*_pfParticles[i0].pt() < 1.0) continue;
        if(pWeight < 0.1) continue;
        PseudoJet curjet( pWeight*_pfParticles[i0].px(), pWeight*_pfParticles[i0].py(), pWeight*_pfParticles[i0].pz(), pWeight*_pfParticles[i0].e());
        curjet.set_user_index(_pfParticles[i0].user_index());            
        particles.push_back(curjet);
    }
    return particles;
}
double puppiContainer::goodVar(PseudoJet &iPart,std::vector<PseudoJet> &iParts, int iOpt) { 
    double Rsub = 0.3; //For tracking using 0.2
                       //double RLrg = 1.0;
    double lPup = 0; 
    lPup = var_within_R(iOpt,iParts,iPart,Rsub);            
    if(iOpt == 6) lPup = lPup * iPart.pt()/pt_within_R(_pfParticles,iPart,Rsub);
    return lPup;
}
void puppiContainer::getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,std::vector<int> &iIsPU,double iQuant,double iPtRMS) { 
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
double puppiContainer::compute(int iOpt,double iVal,double iMed,double iRMS) { 
    if(iOpt == 1 && iVal < iMed) return 0;
    if(iOpt == 1 && iVal > iMed) return 1;
    double lVal = (iVal-iMed)/iRMS;
    //return lVal*fabs(lVal);
    //return TMath::Erf(lVal);//*fabs(lVal);
    return  ROOT::Math::chisquared_cdf(lVal*fabs(lVal),1.);
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
