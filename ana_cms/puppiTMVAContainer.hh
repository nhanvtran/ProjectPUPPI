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
#include "RecoObj.hh"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//......................
class puppiTMVAContainer{
public:
  
  puppiTMVAContainer(std::vector<RecoObj> &inParticles,std::string iWeight="TMVA/weights/TMVAClassificationCategory_PUDisc_v1.weights.xml");  
  ~puppiTMVAContainer();
  std::vector<fastjet::PseudoJet> genParticles(){ return _genParticles; }
  std::vector<fastjet::PseudoJet> pfParticles(){ return _pfParticles; }    
  std::vector<fastjet::PseudoJet> pvParticles(){ return _chargedPV; }        
  std::vector<fastjet::PseudoJet> puParticles(){ return _chargedNoPV; }    
  std::vector<fastjet::PseudoJet> pfchsParticles(){ return _pfchsParticles; }    
  std::vector<fastjet::PseudoJet> puppiEvent     (int iOpt,double iQuant);
  //Load a new event
  void    refresh(std::vector<RecoObj> &inParticles);

protected:
    //Helper Functions
    double  goodVar    (fastjet::PseudoJet &iPart,std::vector<fastjet::PseudoJet> &iParts, int iOpt);    
    void    puppiForAll(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles);
    double  compute  (int iOpt,double iVal,double iMed,double iRMS,double iChi2);
    void    computeMeanRMS();
    std::vector<RecoObj> _recoParticles;
    std::vector<PseudoJet> _pfParticles;
    std::vector<PseudoJet> _pfchsParticles;    
    std::vector<PseudoJet> _genParticles;
    std::vector<PseudoJet> _chargedPV;
    std::vector<PseudoJet> _chargedNoPV;
    std::vector<double> _vals;
    // Global Mean and RMS 
    float fMean;
    float fRMS;
    float fMeanHEta;
    float fRMSHEta;
    // TMVA Variables
    TMVA::Reader *fReader;
    float fTk;
    float fPVtx;
    float fVId;
    float fPt; 
    float fEta; 
    float fDepth;
    float fTime;
    float fPFType;
    float fPuppi;
    float fPuppiO;
    float fPuppiLV;
    float fPuppiOLV;
};

//FASTJET_END_NAMESPACE

