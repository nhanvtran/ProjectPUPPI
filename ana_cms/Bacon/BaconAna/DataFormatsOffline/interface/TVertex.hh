#ifndef BACONANA_DATAFORMATS_TVERTEX_HH
#define BACONANA_DATAFORMATS_TVERTEX_HH

#include <TObject.h>

namespace baconhep
{
  class TVertex : public TObject
  {
    public:
      TVertex():
      nTracksFit(0),
      ndof(0), chi2(0),
      x(0), y(0), z(0)
      {}
      ~TVertex(){}
      
      unsigned int nTracksFit;  // number of tracks in vertex fit
      float        ndof;        // vertex fit number of degrees of freedom
      float        chi2;        // vertex fit chi-square
      float        x,y,z;       // position
    
    ClassDef(TVertex,1)
  };
}
#endif
