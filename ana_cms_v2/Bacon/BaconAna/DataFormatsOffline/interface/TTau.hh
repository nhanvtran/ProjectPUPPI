#ifndef BACONANA_DATAFORMATS_TTAU_HH
#define BACONANA_DATAFORMATS_TTAU_HH

#include "../interface/BaconAnaDefs.hh"
#include <TObject.h>

namespace baconhep
{
  class TTau : public TObject
  {
    public:
      TTau():
      pt(0), eta(0), phi(0), m(0), e(0),
      q(0),
      dzLeadChHad(-999.),
      nSignalChHad(0), nSignalGamma(0),
      ringIso(-1), ringIso2(-1),
      antiEleMVA3(0), antiEleMVA3Cat(0),
      rawIso3Hits(0), rawIsoMVA3(0),
      hpsDisc(0),
      hltMatchBits(0)
      {}
      ~TTau(){}
    
      float         pt, eta, phi, m, e;  // kinematics
      int           q;                   // charge
      float         dzLeadChHad;         // longitudinal impact parameter of leading charged hadron
      unsigned int  nSignalChHad;        // number of charged hadrons in tau
      unsigned int  nSignalGamma;        // number of photons in tau     
      float         ringIso, ringIso2;   // ring isolation
      float         antiEleMVA3;         // anti-electron MVA
      float         antiEleMVA3Cat;      // anti-electron MVA category
      float         rawIso3Hits;         // raw Iso 3 hits
      float         rawIsoMVA3;          // raw Iso MVA3
      unsigned long hpsDisc;             // HPS tau discriminators
      TriggerObjects hltMatchBits;        // HLT matches
         
    ClassDef(TTau,1)
  };
  
  enum EHPSDiscriminator{
    // descriptions from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Discriminators
    kByDecayModeFinding                        = 1UL<<0,   // You will always want to use this (see AN-10-82)
    
    kByVLooseIsolation                         = 1UL<<1,   // isolation cone of 0.3, no PF Charged Candidates with pT > 1.5 GeV/c and no PF Gamma candidates with ET > 2.0 GeV
    kByLooseIsolation                          = 1UL<<2,   // (description N/A)
    kByMediumIsolation                         = 1UL<<3,   // (description N/A)
    kByTightIsolation                          = 1UL<<4,   // (description N/A)
    kByVLooseIsolationDBSumPtCorr              = 1UL<<5,   // (description N/A)
    kByLooseIsolationDBSumPtCorr               = 1UL<<6,   // (description N/A)
    kByMediumIsolationDBSumPtCorr              = 1UL<<7,   // (description N/A)
    kByTightIsolationDBSumPtCorr               = 1UL<<8,   // (description N/A)
    kByVLooseCombinedIsolationDBSumPtCorr      = 1UL<<9,   // isolation cone of 0.3, Delta Beta corrected sum pT of PF charged and PF gamma isolation candidates (pT > 0.5 GeV) less than 3 GeV
    kByLooseCombinedIsolationDBSumPtCorr       = 1UL<<10,  // isolation cone of 0.5, Delta Beta corrected sum pT of PF charged and PF gamma isolation candidates (pT > 0.5 GeV) less than 2 GeV
    kByMediumCombinedIsolationDBSumPtCorr      = 1UL<<11,  // isolation cone of 0.5, Delta Beta corrected sum pT of PF charged and PF gamma isolation candidates (pT > 0.5 GeV) less than 1 GeV
    kByTightCombinedIsolationDBSumPtCorr       = 1UL<<12,  // isolation cone of 0.5, Delta Beta corrected sum pT of PF charged and PF gamma isolation candidates (pT > 0.5 GeV) less than 0.8 GeV
    kByLooseCombinedIsolationDBSumPtCorr3Hits  = 1UL<<13,  // same as ByLooseCombinedIsolationDBSumPtCorr but requiring 3 hits (instead of 8) on track of isolation candidates
    kByMediumCombinedIsolationDBSumPtCorr3Hits = 1UL<<14,  // same as ByMediumCombinedIsolationDBSumPtCorr but requiring 3 hits (instead of 8) on track of isolation candidates
    kByTightCombinedIsolationDBSumPtCorr3Hits  = 1UL<<15,  // same as ByTightCombinedIsolationDBSumPtCorr but requiring 3 hits (instead of 8) on track of isolation candidates
    kByLooseIsolationMVA                       = 1UL<<16,  // BDT based selection using isolation in rings around tau direction and shower shape variables
    kByMediumIsolationMVA                      = 1UL<<17,  //  "
    kByTightIsolationMVA                       = 1UL<<18,  //  "
    kByLooseIsolationMVA2                      = 1UL<<19,  // same as ByLooseIsolationMVA with new training and improved performance
    kByMediumIsolationMVA2                     = 1UL<<20,  // same as ByMediumIsolationMVA with new training and improved performance
    kByTightIsolationMVA2                      = 1UL<<21,  // same as ByTightIsolationMVA with new training and improved performance

    kByLooseElectronRejection                  = 1UL<<22,  // electron pion MVA discriminator < 0.6
    kByMediumElectronRejection                 = 1UL<<23,  // electron pion MVA discriminator < -0.1 and not 1.4442 < |eta| < 1.566
    kByTightElectronRejection                  = 1UL<<24,  // electron pion MVA discriminator < -0.1 and not 1.4442 < |eta| < 1.566 and Brem pattern cuts (see AN-10-387)
    kByMVA3LooseElectronRejection              = 1UL<<25,  // anti-electron MVA discriminator with improved training
    kByMVA3MediumElectronRejection             = 1UL<<26,  //  "
    kByMVA3TightElectronRejection              = 1UL<<27,  //  "
    kByMVA3VTightElectronRejection             = 1UL<<28,  //  "
    
    kByLooseMuonRejection                      = 1UL<<29,  // Tau Lead Track not matched to chamber hits
    kByMediumMuonRejection                     = 1UL<<30,  // Tau Lead Track not matched to global/tracker muon
    kByTightMuonRejection                      = 1UL<<31,  // Tau Lead Track not matched to global/tracker muon and large enough energy deposit in ECAL+HCAL exceeding 0.2 times Lead Track momentum
    kByLooseMuonRejection2                     = 1UL<<32,  // Same as LooseMuonRejection
    kByMediumMuonRejection2                    = 1UL<<33,  // LooseMuonRejection2 && no DT, CSC or RPC Hits in last 2 Stations
    kByTightMuonRejection2                     = 1UL<<34,  // MediumMuonRejection2 && large enough energy deposit in ECAL+HCAL in 1 prong + 0 strip decay mode (SUM(ECAL+HCAL)>0.2*pT)
    kByLooseMuonRejection3                     = 1UL<<35,  // Same as LooseMuonRejection
    kByTightMuonRejection3                     = 1UL<<36   // MediumMuonRejection2 && large enough energy deposit in ECAL+HCAL in 1 prong + 0 strip decay mode (SUM(ECAL+HCAL)>0.2*pT)
  };
}
#endif
