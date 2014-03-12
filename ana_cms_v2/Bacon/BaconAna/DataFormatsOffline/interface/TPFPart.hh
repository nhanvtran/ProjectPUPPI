#ifndef BACONANA_DATAFORMATS_TPFPART_HH
#define BACONANA_DATAFORMATS_TPFPART_HH

#include <TObject.h>


namespace baconhep
{
  class TPFPart : public TObject
  {
    public:
      TPFPart():
	pt(0), eta(0), phi(0), m(0), e(0),
	q(0),pfType(-1),vtxId(-1),
	trkChi2(0),vtxChi2(0),
	ecalE(0),hcalE(0),
	d0(0),dz(0),time(0),depth(0)
    {}
    ~TPFPart(){}
    
    float         pt, eta, phi, m, e;  // kinematics
    int           q;                   // charge
    int           pfType;
    int           vtxId;               // Vertex Id from Vertex Collection
    float         trkChi2;             // Track Chi2
    float         vtxChi2;             // Vertex Chi2
    float         ecalE,hcalE;         // Calorimeter Energy
    float         d0,dz,time,depth;    // Usefule Info
    //float         deltaP;              // Uncertainty on P
    ClassDef(TPFPart,1)
  };
}
#endif
