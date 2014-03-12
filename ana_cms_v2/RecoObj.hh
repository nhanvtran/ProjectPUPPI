#ifndef BACONANA_DATAFORMATS_RECOOBJ_HH
#define BACONANA_DATAFORMATS_RECOOBJ_HH

class RecoObj 
{
public:
      RecoObj():
	pt(0), eta(0), phi(0), m(0),
	pfType(-1),vtxId(-1),
	trkChi2(0),vtxChi2(0),
	id(0),time(0),depth(0)
    {}
    ~RecoObj(){}
    
    float         pt, eta, phi, m;  // kinematics
    int           id;
    int           pfType;
    int           vtxId;               // Vertex Id from Vertex Collection
    float         trkChi2;             // Track Chi2
    float         vtxChi2;             // Vertex Chi2
    float         time,depth;    // Usefule Info
    float         expProb;
    float         expChi2PU;
    float         expChi2;
    float         dZ;
    float         d0;
};
#endif
