#include "NoTrees.hh"
#include "fastjet/internal/base.hh"
#include "fastjet/PseudoJet.hh"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//......................
class puppiContainer{
public:
    // default ctor
    puppiContainer(std::vector<PseudoJet> inParticles, std::vector<int> isPU, std::vector<int> isCh); 
    ~puppiContainer(); 
    
    std::vector<fastjet::PseudoJet> genParticles(){ return _genParticles; }
    std::vector<fastjet::PseudoJet> pfParticles();
    std::vector<fastjet::PseudoJet> pfchsParticles();

    std::vector<float> getPuppiWeights(){ return puppiWeights; };  
    std::vector<float> getPuppiWeights_chLV(){ return puppiWeights_chLV; };  
    std::vector<float> getPuppiWeights_combined(){ return puppiWeights_combined; };  

    std::vector<float> getCleansedWeights(){ return cleansedWeights; };      
    
    std::vector<fastjet::PseudoJet> trimEvent(double vRjet=0.7, double vptcut = 25, double vRsub=0.3, double vfcut=0.05 );
    std::vector<fastjet::PseudoJet> cleanseEvent( double Rsub=0.3 );    
    std::vector<fastjet::PseudoJet> puppiEvent_V1( double Rsub=0.3, double exponent=2. );

    std::vector<fastjet::PseudoJet> puppiEvent( int nPU );

    
protected:
        
    double getAverage(const vector<double> & particles);
    double getRMS(const vector<double> & particles);
            
    double pt_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R);
    double ktWeight_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R, double exponent);    
    
    std::vector<PseudoJet> _pfParticles;
    std::vector<PseudoJet> _pfchsParticles;    
    std::vector<PseudoJet> _genParticles;
    
    std::vector<PseudoJet> _neutrals;
    std::vector<PseudoJet> _chargedLV;    
    std::vector<PseudoJet> _chargedPU;        
    
    std::vector<int> _isPU;
    std::vector<int> _isCh;    
    std::vector<int> _isPFCHS; 
    
    std::vector<float> puppiWeights;
    std::vector<float> puppiWeights_chLV;
    std::vector<float> puppiWeights_combined;
            
    std::vector<float> cleansedWeights;
    
};

//FASTJET_END_NAMESPACE

