#ifndef BACONANA_DATAFORMATS_TEVENTINFO_HH
#define BACONANA_DATAFORMATS_TEVENTINFO_HH

#include "../interface/BaconAnaDefs.hh"
#include <TObject.h>

namespace baconhep 
{
  class TEventInfo : public TObject
  {
    public:
      TEventInfo():
      runNum(0), evtNum(0), lumiSec(0),
      metFilterFailBits(0),
      nPU(0), nPUm(0), nPUp(0),
      nPUmean(0), nPUmeanm(0), nPUmeanp(0),
      pvx(0), pvy(0), pvz(0),
      bsx(0), bsy(0), bsz(0), 
      pfMET(0), pfMETphi(0), pfMETCov00(0), pfMETCov01(0), pfMETCov11(0),
      mvaMET(0), mvaMETphi(0), mvaMETCov00(0), mvaMETCov01(0), mvaMETCov11(0),
      mvaMETU(0), mvaMETUphi(0), mvaMETUCov00(0), mvaMETUCov01(0), mvaMETUCov11(0),
      trkMET(0), trkMETphi(0),
      rhoIso(0), rhoJet(0),
      triggerBits(0),
      hasGoodPV(false)             
      {}
      ~TEventInfo(){}

      unsigned int  runNum, evtNum, lumiSec;                   // run number, event number, lumi section in data
      unsigned int  metFilterFailBits;                         // MET filter fail bits
      unsigned int  nPU, nPUm, nPUp;                           // number of reconstructed pile-up vertices
      float         nPUmean, nPUmeanm, nPUmeanp;               // expected number of reconstructed pile-up vertices
      float         pvx, pvy, pvz;                             // best primary vertex
      float         bsx, bsy, bsz;                             // beamspot
      float         pfMET, pfMETphi;                           // particle flow MET
      float         pfMETCov00, pfMETCov01, pfMETCov11;
      float         mvaMET, mvaMETphi;                         // MVA MET
      float         mvaMETCov00, mvaMETCov01, mvaMETCov11;
      float         mvaMETU, mvaMETUphi;
      float         mvaMETUCov00, mvaMETUCov01, mvaMETUCov11;
      float         trkMET, trkMETphi;                         // track MET
      float         rhoIso, rhoJet;                            // event energy density
      TriggerBits   triggerBits;                               // fired trigger bits
      bool          hasGoodPV;                                 // event has a good PV?

    ClassDef(TEventInfo,1)
  };
}
#endif
