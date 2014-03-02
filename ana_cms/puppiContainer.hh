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
  puppiContainer(std::vector<PseudoJet> inParticles, std::vector<int> isPU, std::vector<int> isCh,bool iDiscretize); 
    ~puppiContainer(); 
    
    std::vector<fastjet::PseudoJet> genParticles(){ return _genParticles; }
    std::vector<fastjet::PseudoJet> pfParticles(){ return _pfParticles; }    
    std::vector<fastjet::PseudoJet> pvParticles(){ return _chargedPV; }        
    std::vector<fastjet::PseudoJet> puParticles(){ return _chargedNoPV; }    
    std::vector<fastjet::PseudoJet> pfchsParticles(){ return _pfchsParticles; }    
    void discretize(std::vector<fastjet::PseudoJet> &discreteParticles,std::vector<fastjet::PseudoJet> &iParticles);
    double goodVar(fastjet::PseudoJet &iPart,std::vector<fastjet::PseudoJet> &iParts, int iOpt);    

protected:
        
    std::vector<fastjet::PseudoJet> trimEvent();
    std::vector<PseudoJet> _pfParticles;
    std::vector<PseudoJet> _pfchsParticles;    
    std::vector<PseudoJet> _genParticles;
    std::vector<PseudoJet> _chargedPV;
    std::vector<PseudoJet> _chargedNoPV;
    std::vector<int> _isPU;
    std::vector<int> _isCh;    
    
};

//FASTJET_END_NAMESPACE

