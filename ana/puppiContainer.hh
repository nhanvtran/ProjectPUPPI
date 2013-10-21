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
    std::vector<fastjet::PseudoJet> pfParticles(){ return _pfParticles; }    
    std::vector<fastjet::PseudoJet> pfchsParticles(){ return _pfchsParticles; }  
    std::vector<float> getPuppiWeights(){ return puppiWeights; };  
    std::vector<float> getCleansedWeights(){ return cleansedWeights; };      
    
protected:
        
    std::vector<fastjet::PseudoJet> trimEvent(double vRjet=0.7, double vptcut = 25, double vRsub=0.3, double vfcut=0.05 );
    std::vector<fastjet::PseudoJet> cleanseEvent( double Rsub=0.3 );    
    std::vector<fastjet::PseudoJet> puppiEvent( double Rsub=0.3, double exponent=2. );

    double pt_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R);
    double ktWeight_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R, double exponent);    
    
    std::vector<PseudoJet> _pfParticles;
    std::vector<PseudoJet> _pfchsParticles;    
    std::vector<PseudoJet> _genParticles;
    std::vector<int> _isPU;
    std::vector<int> _isCh;    
    
    std::vector<float> puppiWeights;
    std::vector<float> cleansedWeights;
    
};

//FASTJET_END_NAMESPACE

