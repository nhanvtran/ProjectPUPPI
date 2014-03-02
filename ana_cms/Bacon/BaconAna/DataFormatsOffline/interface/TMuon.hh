#ifndef BACONANA_DATAFORMATS_TMUON_HH
#define BACONANA_DATAFORMATS_TMUON_HH

#include "../interface/BaconAnaDefs.hh"
#include <TObject.h>

namespace baconhep
{
  class TMuon : public TObject
  {
    public:
      TMuon():
      pt(0), eta(0), phi(0), ptErr(0), ptHZZ4l(0),
      staPt(0), staEta(0), staPhi(0),
      pfPt(0), pfEta(0), pfPhi(0),
      trkIso03(-1), ecalIso03(-1), hcalIso03(-1),
      chHadIso03(-1), gammaIso03(-1), neuHadIso03(-1), puIso03(-1),
      chHadIso04(-1), gammaIso04(-1), neuHadIso04(-1), puIso04(-1),
      d0(-999.), dz(-999.), sip3d(-999.),
      tkNchi2(-999.), muNchi2(-999.),
      trkKink(0), glbKink(0),
      q(0),
      nValidHits(0),
      typeBits(0), selectorBits(0),
      nTkHits(0), nPixHits(0),
      nTkLayers(0), nPixLayers(0),
      nMatchStn(0),
      trkID(-1),
      hltMatchBits(0)
      {}
      ~TMuon(){}
    
      float         pt, eta, phi, ptErr;                           // kinematics
      float         ptHZZ4l;                                       // HZZ4l lepton corrections
      float         staPt, staEta, staPhi;                         // STA track kinematics
      float         pfPt, pfEta, pfPhi;                            // matched PFCandidate
      float         trkIso03, ecalIso03, hcalIso03;                // detector isolation
      float         chHadIso03, gammaIso03, neuHadIso03, puIso03;  // PF isolation variables
      float         chHadIso04, gammaIso04, neuHadIso04, puIso04;
      float         d0, dz, sip3d;                                 // impact parameter
      float         tkNchi2, muNchi2;                              // track fit normalized chi-square
      float         trkKink, glbKink;                              // track kink
      int           q;                                             // charge
      int           nValidHits;                                    // number of valid muon hits in global fit
      unsigned int  typeBits;                                      // muon type bits
      unsigned int  selectorBits;                                  // MuonSelector bits
      unsigned int  nTkHits, nPixHits;                             // number of hits in tracker
      unsigned int  nTkLayers, nPixLayers;                         // number of hit layers in tracker
      unsigned int  nMatchStn;                                     // number of stations with muon segments
      int           trkID;                                         // tracker track ID (unique per event)
      TriggerObjects hltMatchBits;                                  // HLT matching
          
    ClassDef(TMuon,1)
  };

  enum EMuType
  {
    kGlobal     = 2,
    kTracker    = 4,
    kStandalone = 8,
    kCaloMuon   = 16,
    kPFMuon     = 32,
    kRPCMuon    = 64
  };

  enum EMuSelectorBit
  { 
    // descriptions from DataFormats/MuonReco/interface/MuonSelectors.h
    kAll  			      = 0x0000001,  // dummy options - always true
    kAllGlobalMuons		      = 0x0000002,  // checks isGlobalMuon flag
    kAllStandAloneMuons		      = 0x0000004,  // checks isStandAloneMuon flag
    kAllTrackerMuons		      = 0x0000008,  // checks isTrackerMuon flag
    kTrackerMuonArbitrated	      = 0x0000010,  // resolve ambiguity of sharing segments
    kAllArbitrated		      = 0x0000020,  // all muons with the tracker muon arbitrated
    kGlobalMuonPromptTight	      = 0x0000040,  // global muons with tighter fit requirements
    kTMLastStationLoose		      = 0x0000080,  // penetration depth loose selector
    kTMLastStationTight		      = 0x0000100,  // penetration depth tight selector
    kTM2DCompatibilityLoose	      = 0x0000200,  // likelihood based loose selector
    kTM2DCompatibilityTight	      = 0x0000400,  // likelihood based tight selector
    kTMOneStationLoose		      = 0x0000800,  // require one well matched segment
    kTMOneStationTight		      = 0x0001000,  // require one well matched segment
    kTMLastStationOptimizedLowPtLoose = 0x0002000,  // combination of TMLastStation and TMOneStation
    kTMLastStationOptimizedLowPtTight = 0x0004000,  // combination of TMLastStation and TMOneStation
    kGMTkChiCompatibility 	      = 0x0008000,  // require tk stub have good chi2 relative to glb track
    kGMStaChiCompatibility	      = 0x0010000,  // require sta stub have good chi2 compatibility relative to glb track
    kGMTkKinkTight		      = 0x0020000,  // require a small kink value in the tracker stub
    kTMLastStationAngLoose	      = 0x0040000,  // TMLastStationLoose with additional angular cuts
    kTMLastStationAngTight	      = 0x0080000,  // TMLastStationTight with additional angular cuts
    kTMOneStationAngLoose 	      = 0x0100000,  // TMOneStationLoose with additional angular cuts
    kTMOneStationAngTight 	      = 0x0200000,  // TMOneStationTight with additional angular cuts
    //The two algorithms that follow are identical to what were known as
    //TMLastStationOptimizedLowPt* (sans the Barrel) as late as revision
    //1.7 of this file. The names were changed because indeed the low pt
    //optimization applies only to the barrel region, whereas the sel-
    //ectors above are more efficient at low pt in the endcaps, which is
    //what we feel is more suggestive of the algorithm name. This will be
    //less confusing for future generations of CMS members, I hope...
    kTMLastStationOptimizedBarrelLowPtLoose = 0x0400000,  // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
    kTMLastStationOptimizedBarrelLowPtTight = 0x0800000,  // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
    kRPCMuLoose                             = 0x1000000   // checks isRPCMuon flag (require two well matched hits in different RPC layers)
  };
}
#endif
