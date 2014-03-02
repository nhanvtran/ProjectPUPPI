#include "Bacon/BaconAna/DataFormatsOffline/interface/TPFPart.hh"

#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TLorentzVector.h"

class GenMatcher { 
public:
  GenMatcher(TTree *iTree);
  ~GenMatcher();
  void setup(TTree *iTree);
  bool load(int iEvent);
  void match();
  void fillGen();

protected: 
  double deltaR(float iEta0,float iPhi0,float iEta1,float iPhi1);
  float  matchDilution(int iId,float iPt, float iEta,float iPhi);
  TLorentzVector matchGen(float iEta,float iPhi,float iDilution);
  void   fillVars(TLorentzVector iMatch,baconhep::TPFPart *iPart);
  void   reset();
  TClonesArray *fPFPart;
  TClonesArray *fGenPart;
  TBranch      *fPFPartBr;
  TBranch      *fGenPartBr;
  TTree        *fTree;
  float fPtGen;
  float fEtaGen;
  float fPhiGen;
  float fMGen;
  float fDRGen;
  float fPtRec;
  float fEtaRec;
  float fPhiRec;
  float fMRec;
  float fTrk;
  float fVtx;
  int   fVtxId;
};
