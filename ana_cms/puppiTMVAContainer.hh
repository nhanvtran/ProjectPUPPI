#include "NoTrees.hh"
#include "fastjet/internal/base.hh"
#include "fastjet/PseudoJet.hh"
#include "TSystem.h"
#include "TROOT.h"
//#include "TMVA/TMVA/TMVAGui.C"
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
//#include "TMVA/MethodCategory.h"
#include "TMVA/Tools.h"
#endif


using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//......................
class puppiTMVAContainer{
public:
    // default ctor
  puppiTMVAContainer(std::vector<PseudoJet> &inParticles, std::vector<int> &isPU, std::vector<int> &isCh, std::vector<float> &iTk, std::vector<float> &iVtx,std::vector<int> &iVId,std::string iWeight="TMVA/weights/TMVAClassificationCategory_PUDisc_v1.weights.xml",bool iVtx=true,bool iDiscretize=false,bool iCHS=false); 
    ~puppiTMVAContainer(); 
  void discretize(std::vector<fastjet::PseudoJet> &discreteParticles,std::vector<fastjet::PseudoJet> &iParticles,bool iPtCut=false);
    std::vector<fastjet::PseudoJet> genParticles(){ return _genParticles; }
    std::vector<fastjet::PseudoJet> pfParticles(){ return _pfParticles; }    
    std::vector<fastjet::PseudoJet> pvParticles(){ return _chargedPV; }        
    std::vector<fastjet::PseudoJet> puParticles(){ return _chargedNoPV; }    
    std::vector<fastjet::PseudoJet> pfchsParticles(){ return _pfchsParticles; }    
    std::vector<fastjet::PseudoJet> puppiEvent     (int iOpt,double iQuant);
  //Load a new event
  void    refresh(std::vector<fastjet::PseudoJet> &inParticles, std::vector<int> &isPU, std::vector<int> &isCh,std::vector<float> &iTk, std::vector<float> &iVtx, std::vector<int> &iVId);

protected:
        
    //Helper Functions
    double  goodVar  (fastjet::PseudoJet &iPart,std::vector<fastjet::PseudoJet> &iParts, int iOpt);    
    void    getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,std::vector<int> &iIsPU,double iQuant);
    double  compute  (int iOpt,double iVal,double iMed,double iRMS);
    std::vector<PseudoJet> _pfParticles;
    std::vector<PseudoJet> _pfchsParticles;    
    std::vector<PseudoJet> _genParticles;
    std::vector<PseudoJet> _chargedPV;
    std::vector<PseudoJet> _chargedNoPV;
    std::vector<double> _vals;
    std::vector<int>   _isPU;
    std::vector<int>   _isCh;    
    std::vector<float> _vtx;
    std::vector<float> _tk;    
    std::vector<int>   _vid;    
    double fMed;
    double fRMS;
    bool  fCHS;
    bool  fVtx;
    bool  fDiscretize;
    TMVA::Reader *fReader;
    TMVA::Reader *fReaderNoVtx;
    float fTk;
    float fPVtx;
    float fVId;
    float fPt; 
    float fEta; 
    float fPuppi;
    float fPuppiO;
    float fPuppiLV;
    float fPuppiOLV;
};

//FASTJET_END_NAMESPACE

