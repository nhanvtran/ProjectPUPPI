#ifndef BACONANA_DATAFORMATS_TJET_HH
#define BACONANA_DATAFORMATS_TJET_HH

#include "../interface/BaconAnaDefs.hh"
#include <TObject.h>

namespace baconhep
{
  class TJet : public TObject
  {
    public:
      TJet():
      pt(0), eta(0), phi(0), mass(0), ptRaw(0),unc(0),
      area(0), d0(-999),dz(-999),
      csv(-2),csv1(-2),csv2(-2),mva(-2),
      qgid(-2),qg1(-2),qg2(-2),
      tau1(-1),tau2(-1),tau3(-1),tau4(-1),
      prunedm(0),
      nCharged(0), nNeutrals(0), nParticles(0),
      beta(0),betaStar(0),dR2Mean(0),ptD(0),q(0),pull(0),pullAngle(0),
      chEmFrac(0), neuEmFrac(0), chHadFrac(0), neuHadFrac(0),
      mcFlavor(0), 
      genpt(0),geneta(0),genphi(0),genm(0),
      hltMatchBits(0)
      {}
      ~TJet(){}
    
      float pt, eta, phi, mass, ptRaw,unc;                   // kinematics
      float area;                                        // jet area (from FastJet)
      float d0, dz;                                      // impact parameter of leading charged constituent
      float csv,csv1,csv2;                               // combined secondary vertex based b-tagger and for two leading subjets
      float mva;                                         // Pile-up discriminator MVA
      float qgid,qg1,qg2;                                // Quark Gluon Id for jet and subjets
      float tau1,tau2,tau3,tau4;                         // N Subjettiness variables
      float prunedm;                                     // Pruned Jet Mass
      unsigned int nCharged, nNeutrals, nParticles;      // constituent multiplicity
      float beta,betaStar,dR2Mean;                       // PU Jet Id Variables
      float ptD,q,pull,pullAngle;                        // Additional color flow variables
      float chEmFrac, neuEmFrac, chHadFrac, neuHadFrac;  // fractional contribution by type
      int mcFlavor,mcFlavorPhys;                                      //Flavor
      float genpt,geneta,genphi,genm;                    // Gen Matched values
      TriggerObjects hltMatchBits;                        // HLT matching
          
    ClassDef(TJet,1)
  };
}
#endif
