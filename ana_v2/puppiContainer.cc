#include "puppiContainer.hh"

#include "fastjet/internal/base.hh"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


puppiContainer::puppiContainer(std::vector<fastjet::PseudoJet> inParticles, std::vector<int> isPU, std::vector<int> isCh){
    
    _pfParticles = inParticles;
    _isPU = isPU;
    _isCh = isCh;
    
    _genParticles.resize(0);
    _pfchsParticles.resize(0);    
    _neutrals.resize(0);
    _chargedLV.resize(0);    
    _chargedPU.resize(0);    
        
    double etaTracker = 5.;
    
    for (unsigned int i = 0; i < _pfParticles.size(); i++){
        float weightGEN = 0.;
        float weightPFCHS = 0.;
    
        // fill vector of pseudojets
        if (fabs(_pfParticles[i].eta()) < 5){        
            
            PseudoJet curjet(_pfParticles[i].px(), _pfParticles[i].py(), _pfParticles[i].pz(), _pfParticles[i].e());
            
            if (isPU[i] == 0){
                weightGEN = 1.;
            }
            
            if ((isPU[i] == 0) || (isPU[i] == 1 && isCh[i] == 0 && fabs(_pfParticles[i].eta()) < etaTracker) || (isPU[i] == 1 && fabs(_pfParticles[i].eta()) > etaTracker)){
                weightPFCHS = 1.;
            }
            if (weightGEN > 0){
                PseudoJet curjetGEN( weightGEN*_pfParticles[i].px(), weightGEN*_pfParticles[i].py(), weightGEN*_pfParticles[i].pz(), weightGEN*_pfParticles[i].e());
                _genParticles.push_back( curjetGEN );
            }
            if (weightPFCHS > 0){
                PseudoJet curjetPFCHS( weightPFCHS*_pfParticles[i].px(), weightPFCHS*_pfParticles[i].py(), weightPFCHS*_pfParticles[i].pz(), weightPFCHS*_pfParticles[i].e());            
                _pfchsParticles.push_back( curjetPFCHS );
                
                if (isCh[i] == 0){
                    _neutrals.push_back( curjet );
                }
                else{
                    _chargedLV.push_back( curjet );
                }
            }
            if (weightPFCHS == 0){
                _chargedPU.push_back( curjet );
            }
        }
    }
    
    puppiWeights.resize(0);
    cleansedWeights.resize(0);    
}

puppiContainer::~puppiContainer(){}

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
    puppiWeights.resize(0);
    
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
        puppiWeights.push_back( weight1 );    
    }
    
    return answer;
}

std::vector<fastjet::PseudoJet> puppiContainer::puppiEvent(int iPU){
    
    double Rsub = 0.3;///Using 0.3 radius
    std::vector<PseudoJet> answer;
    answer.resize(0);
    // for shape puppi
    std::vector<double> lVals;
    double lAvg  = 0; 
    // for charged LV puppi
    std::vector<double> lVals_chLV;
    double lAvg_chLV  = 0; 
    // for charged PU puppi
    std::vector<double> lVals_chPU;
    double lAvg_chPU  = 0; 
    
    for(unsigned int i=0; i<_pfParticles.size(); i++){
        //if((!_isPU[i] && iPU==0) || ( _isPU[i] && iPU==1)) {
        //    lVals.push_back(0.);
        //    continue;
        //}
        double puppi_Rsub = puppi_within_R(_pfParticles,_pfParticles[i],Rsub,false);            
        double puppi_Rsub_chLV = puppi_within_R(_chargedLV,_pfParticles[i],Rsub,false);            
        double puppi_Rsub_chPU = puppi_within_R(_chargedPU,_pfParticles[i],Rsub,false);            

        //std::cout << "puppi_Rsub = " << puppi_Rsub << std::endl;
        
        lVals.push_back(puppi_Rsub);
        lAvg += lVals[i];            

        lVals_chLV.push_back(puppi_Rsub_chLV);
        lAvg_chLV += lVals_chLV[i];            

        lVals_chPU.push_back(puppi_Rsub_chPU);
        lAvg_chPU += lVals_chPU[i];            
        
        //Was applying a pt cut
        //if(_pfParticles[i].pt() > 0.2) { 
        //lVals.push_back(puppi_Rsub);
        //    lAvg  += puppi_Rsub;
        // }
    }
    
    std::cout << "lVals.size() = " << lVals.size() << ", _pfParticles.size() = " << _pfParticles.size() << ", lAvg = " << lAvg << std::endl;
    lAvg /= lVals.size();
    double lRMS = 0; 
    for(int i0 = 0 ;i0 < lVals.size(); i0++) {lRMS += (lVals[i0]-lAvg)*(lVals[i0]-lAvg);}
    lRMS/=lVals.size(); 
    std::cout << "Avg = " << lAvg << ", RMS = " << lRMS << std::endl;

    lAvg_chLV /= lVals_chLV.size();
    double lRMS_chLV = 0; 
    for(int i0 = 0 ;i0 < lVals_chLV.size(); i0++) {lRMS_chLV += (lVals_chLV[i0]-lAvg_chLV)*(lVals_chLV[i0]-lAvg_chLV);}
    lRMS_chLV/=lVals_chLV.size(); 

    lAvg_chPU /= lVals_chPU.size();
    double lRMS_chPU = 0; 
    for(int i0 = 0 ;i0 < lVals_chPU.size(); i0++) {lRMS_chPU += (lVals_chPU[i0]-lAvg_chPU)*(lVals_chPU[i0]-lAvg_chPU);}
    lRMS_chPU/=lVals_chPU.size(); 
    
    for(int i0 = 0; i0 < lVals.size(); i0++) { 
        double pWeight = lVals[i0]; 
        double pWeight_chLV = lVals_chLV[i0]; 
        double pWeight_chPU = lVals_chPU[i0]; 
        
        pWeight = erf(-(lVals[i0] - lAvg) /sqrt(lRMS));
        pWeight_chLV = erf(-(lVals_chLV[i0] - lAvg_chLV) /sqrt(lRMS_chLV));
        pWeight_chPU = erf(-(lVals_chPU[i0] - lAvg_chPU) /sqrt(lRMS_chPU));
        
        //std::cout << "pWeight = " << pWeight << std::endl;
        PseudoJet curjet( pWeight*_pfParticles[i0].px(), pWeight*_pfParticles[i0].py(), pWeight*_pfParticles[i0].pz(), pWeight*_pfParticles[i0].e());
        answer.push_back( curjet );
        puppiWeights.push_back( pWeight );
        puppiWeights_chLV.push_back( pWeight_chLV );
        puppiWeights_chPU.push_back( pWeight_chPU );
    }
    return answer; 

}


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
//FASTJET_END_NAMESPACE

