#include "NoTrees.hh"
#include "fastjet/internal/base.hh"
#include "fastjet/PseudoJet.hh"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//......................
class puppiCleanContainer{
public:
    // default ctor
    puppiCleanContainer(std::vector<PseudoJet> inParticles, std::vector<int> isPU, std::vector<int> isCh); 
    ~puppiCleanContainer(); 
    
    std::vector<fastjet::PseudoJet> genParticles(){ return _genParticles; }
    std::vector<fastjet::PseudoJet> pfParticles(){ return _pfParticles; }    
    std::vector<fastjet::PseudoJet> pvParticles(){ return _chargedPV; }        
    std::vector<fastjet::PseudoJet> puParticles(){ return _chargedNoPV; }    
    std::vector<fastjet::PseudoJet> pfchsParticles(){ return _pfchsParticles; }    
    
protected:
        
    std::vector<fastjet::PseudoJet> trimEvent();
    std::vector<fastjet::PseudoJet> puppiEvent     (int iOpt,double iQuant);
    //Helper Functions
    double  goodVar  (fastjet::PseudoJet &iPart,std::vector<fastjet::PseudoJet> &iParts, int iOpt);    
  void   getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,std::vector<int> &iIsPU,double iQuant);
    double  compute  (int iOpt,double iVal,double iMed,double iRMS);
    std::vector<PseudoJet> _pfParticles;
    std::vector<PseudoJet> _pfchsParticles;    
    std::vector<PseudoJet> _genParticles;
    std::vector<PseudoJet> _chargedPV;
    std::vector<PseudoJet> _chargedNoPV;
    std::vector<double> _vals;
    std::vector<int> _isPU;
    std::vector<int> _isCh;    
    double fMed;
    double fRMS;
};

//FASTJET_END_NAMESPACE

