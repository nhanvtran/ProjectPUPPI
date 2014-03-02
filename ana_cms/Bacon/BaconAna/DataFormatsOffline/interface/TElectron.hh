#ifndef BACONANA_DATAFORMATS_TELECTRON_HH
#define BACONANA_DATAFORMATS_TELECTRON_HH

#include "../interface/BaconAnaDefs.hh"
#include <TObject.h>

namespace baconhep
{
  class TElectron : public TObject
  {
    public:
      TElectron():
      pt(0), eta(0), phi(0),
      scEt(0), scEta(0), scPhi(0),
      ptHZZ4l(0), ptErrHZZ4l(0), scEtHZZ4l(0),
      r9(0), ecalEnergy(0),
      pfPt(0), pfEta(0), pfPhi(0),
      trkIso03(-1), ecalIso03(-1), hcalIso03(-1),
      chHadIso03(-1), gammaIso03(-1), neuHadIso03(-1),      
      chHadIso04(-1), gammaIso04(-1), neuHadIso04(-1),      
      d0(-999.), dz(-999.), sip3d(-999.),
      sieie(0), eoverp(0), hovere(0), fbrem(0),
      dEtaIn(0), dPhiIn(0),
      mva(-999.),
      q(0), classification(-999.),
      isConv(false),
      nMissingHits(0),
      typeBits(0), fiducialBits(0),
      scID(-1), trkID(-1),
      hltMatchBits(0)
      {}
      ~TElectron(){}
    
      float         pt, eta, phi;                         // kinematics
      float         scEt, scEta, scPhi;                   // supercluster kinematics
      float         ptHZZ4l, ptErrHZZ4l, scEtHZZ4l;       // HZZ4l lepton corrections
      float         r9, ecalEnergy;                       // ECAL energy quantities
      float         pfPt, pfEta, pfPhi;                   // matching PF-candidate kinematics
      float         trkIso03, ecalIso03, hcalIso03;       // detector isolation
      float         chHadIso03, gammaIso03, neuHadIso03;  // PF isolation variables
      float         chHadIso04, gammaIso04, neuHadIso04;
      float         d0, dz, sip3d;                        // impact parameter
      float         sieie;                                // shower shape
      float         eoverp;                               // E/p
      float         hovere;                               // H/E
      float         fbrem;                                // brem fraction
      float         dEtaIn, dPhiIn;                       // track-supercluster matching
      float         mva;                                  // electron ID MVA value
      int           q;                                    // charge
      int           classification;                       // electron classification
      bool          isConv;                               // identified by track fit based conversion finder?
      unsigned int  nMissingHits;                         // number of missing expected inner hits 
      unsigned int  typeBits;                             // electron type
      unsigned int  fiducialBits;                         // ECAL fiducial region bits
      int           scID;                                 // supercluster ID number (unique per event)
      int           trkID;                                // track ID number (unique per event)
      TriggerObjects hltMatchBits;                         // HLT matches
      
    ClassDef(TElectron,1)
  };

  enum EEleType
  {
    kEcalDriven    = 1,
    kTrackerDriven = 2
  };
}
#endif
