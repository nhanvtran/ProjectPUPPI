#include "GenMatcher.hh"
#include "Bacon/BaconAna/DataFormatsOffline/interface/TGenParticle.hh"

#include "TFile.h"
#include "TMath.h"
#include <iostream>

using namespace baconhep;

GenMatcher::GenMatcher(TTree *iTree) { 
  fPFPart  = new TClonesArray("baconhep::TPFPart");
  fGenPart = new TClonesArray("baconhep::TGenParticle");
  iTree->SetBranchAddress("PFPart",       &fPFPart);
  iTree->SetBranchAddress("GenParticle",  &fGenPart);
  fPFPartBr  = iTree->GetBranch("PFPart");
  fGenPartBr = iTree->GetBranch("GenParticle");
}
GenMatcher::~GenMatcher()  { 
  delete fPFPart;
  delete fGenPart;
  delete fPFPartBr;
  delete fGenPartBr;
}
void GenMatcher::setup(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("genpt" ,&fPtGen ,"fPtGen/F");
  fTree->Branch("geneta",&fEtaGen,"fEtaGen/F");
  fTree->Branch("genhi" ,&fPhiGen,"fPhiGen/F");
  fTree->Branch("genm"  ,&fMGen  ,"fMGen/F");
  fTree->Branch("gendr" ,&fDRGen ,"fDRGen/F");
  fTree->Branch("pt"    ,&fPtRec ,"fPtRec/F");
  fTree->Branch("eta"   ,&fEtaRec,"fEtaRec/F");
  fTree->Branch("phi"   ,&fPhiRec,"fPhiRec/F");
  fTree->Branch("m"     ,&fMRec  ,"fMRec/F");
  fTree->Branch("trk"   ,&fTrk   ,"fTrk/F");
  fTree->Branch("vtx"   ,&fVtx   ,"fVtx/F");
  fTree->Branch("vtxid" ,&fVtxId ,"fVtxId/I");
  fTree->Branch("time"  ,&fTime  ,"fTime/F");
  fTree->Branch("depth" ,&fDepth ,"fDepth/F");
  fTree->Branch("pftype",&fPFType,"fPFType/I");
  fTree->Branch("d0"    ,&fD0    ,"fD0/F");
  fTree->Branch("dZ"    ,&fDZ    ,"fDZ/F");
}
bool GenMatcher::load(int iEvent) { 
  fPFPart   ->Clear();
  fGenPart  ->Clear();
  fPFPartBr ->GetEntry(iEvent);
  fGenPartBr->GetEntry(iEvent);
}
float GenMatcher::matchDilution(int iId,float iPt, float iEta,float iPhi) { 
  float lPtTot = iPt;
  for(int i0 = 0; i0 < fPFPart->GetEntriesFast(); i0++) { 
    if(i0 == iId) continue;
    TPFPart *pPart = (TPFPart*)((*fPFPart)[i0]);
    double pDR = deltaR(iEta,iPhi,pPart->eta,pPart->phi);
    if(pDR > 0.05+0.08*(fabs(iEta) > 2.5)) continue;
    lPtTot += pPart->pt;
  }
  return iPt/lPtTot;
}
TLorentzVector GenMatcher::matchGen(float iEta,float iPhi,float iDilution) { 
  TLorentzVector lVec;
  for(int i0 = 0; i0 < fGenPart->GetEntriesFast(); i0++) { 
    TGenParticle *lPart = (TGenParticle*)((*fGenPart)[i0]);
    if(lPart->status != 1) continue;
    if(fabs(lPart->pdgId) == 12 || 
       fabs(lPart->pdgId) == 14 || 
       fabs(lPart->pdgId) == 16) continue;
    double pDR = deltaR(iEta,iPhi,lPart->eta,lPart->phi);
    if(pDR > 0.05 + 0.08*(fabs(iEta) > 2.5)) continue;
    TLorentzVector pVec; pVec.SetPtEtaPhiM(lPart->pt*iDilution,lPart->eta,lPart->phi,lPart->mass*iDilution); //Dilute particle by energy
    lVec += pVec;
  }
  return lVec;
} 
void GenMatcher::match() { 
  reset(); fTree->Fill();
  for  (int i0 = 0; i0 < fPFPart->GetEntriesFast(); i0++) { 
    TPFPart *pPart = (TPFPart*)((*fPFPart)[i0]);
    float pDilution       = matchDilution(i0,pPart->pt,pPart->eta,pPart->phi);
    TLorentzVector pMatch = matchGen(pPart->eta,pPart->phi,pDilution);
    fillVars(pMatch,pPart);
  }
}
void GenMatcher::fillGen() { 
  reset(); fTree->Fill();
  for  (int i0 = 0; i0 < fGenPart->GetEntriesFast(); i0++) { 
    TGenParticle *pPart = (TGenParticle*)((*fGenPart)[i0]);
    if(pPart->status != 1) continue; 
    if(fabs(pPart->pdgId) == 12 || 
       fabs(pPart->pdgId) == 14 || 
       fabs(pPart->pdgId) == 16) continue;
    fPtGen  = pPart->pt;
    fEtaGen = pPart->eta;
    fPhiGen = pPart->phi;
    fMGen   = pPart->mass;
    if(fPtGen < 0) continue;
    fTree->Fill();
  }
}
void GenMatcher::fillVars(TLorentzVector iMatch,TPFPart *iPart) { 
  reset();
  if(iMatch.Pt() > 0) { 
    fPtGen  = iMatch.Pt(); 
    fEtaGen = iMatch.Eta();
    fPhiGen = iMatch.Phi(); 
    fMGen   = iMatch.M();
    fDRGen  = deltaR(iPart->eta,iPart->phi,iMatch.Eta(),iMatch.Phi());
  }
  fPtRec  = iPart->pt; 
  fEtaRec = iPart->eta;
  fPhiRec = iPart->phi; 
  fMRec   = iPart->m;
  fTrk    = iPart->trkChi2;
  fVtx    = iPart->vtxChi2;
  fVtxId  = iPart->vtxId;
  fDepth  = iPart->depth;
  fTime   = iPart->time;
  fPFType = iPart->pfType;
  fD0     = iPart->d0;
  fDZ     = iPart->dz;
  fTree->Fill();
}
double GenMatcher::deltaR(float iEta0,float iPhi0,float iEta1,float iPhi1) { 
  float lDPhi = fabs(iPhi0-iPhi1); 
  if(lDPhi > 2.*TMath::Pi()-lDPhi) lDPhi = 2.*TMath::Pi()-lDPhi;
  float lDEta = iEta0-iEta1;
  return sqrt(lDPhi*lDPhi+lDEta*lDEta);
}
void GenMatcher::reset() { 
  fPtGen  = -1; 
  fEtaGen = -1;
  fPhiGen = -1; 
  fMGen   = -1;
  fDRGen  = -1;
  fPtRec  = -1; 
  fEtaRec = -1;
  fPhiRec = -1; 
  fMRec   = -1;
  fTrk    = -1;
  fVtx    = -1;
  fVtxId  = -1;
  fTime   = -1;
  fDepth  = -1;
  fPFType = -1;
  fD0     = -1;
  fDZ     = -1; 
}
