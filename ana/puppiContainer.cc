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
    
    for (unsigned int i = 0; i < _pfParticles.size(); i++){

        // fill vector of pseudojets
        if (fabs(_pfParticles[i].eta()) < 5){        
            if (isPU[i] == 0){
                _genParticles.push_back( _pfParticles[i] );            
            }
            if ((isPU[i] == 0) || (isPU[i] == 1 && isCh[i] == 0 && fabs(_pfParticles[i].eta()) < 2.5) || (isPU[i] == 1 && fabs(_pfParticles[i].eta()) > 2.5)){
                _pfchsParticles.push_back( _pfParticles[i] );
            }
        }

    }
    
}

puppiContainer::~puppiContainer(){}

std::vector<fastjet::PseudoJet> puppiContainer::trimEvent(){
    
    //std::cout << "trimming event..." << std::endl;
    std::vector<PseudoJet> answer;
    answer.resize(0);
    
    // -- event trimming parameters -- 
    double Rjet = 1;
    double ptcut = 25;
    double Rsub = 0.3;
    double fcut = 0.05;	
    
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

std::vector<fastjet::PseudoJet> puppiContainer::cleanseEvent(){

    //std::cout << "cleansing event..." << std::endl;
    std::vector<PseudoJet> answer;
    answer.resize(0);

    // -- event cleansing parameters -- 
    double Rsub_cleansing = 0.3;

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
        else weight1 = 0.;

        if (weight1 != 0.0){
            PseudoJet curjet( weight1*_pfParticles[i].px(), weight1*_pfParticles[i].py(), weight1*_pfParticles[i].pz(), weight1*_pfParticles[i].e());
            answer.push_back( curjet );
        }
    }
    
    return answer;
}

//FASTJET_END_NAMESPACE

