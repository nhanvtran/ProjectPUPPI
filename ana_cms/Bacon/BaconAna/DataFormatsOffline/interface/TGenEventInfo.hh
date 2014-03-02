#ifndef BACONANA_DATAFORMATS_TGENEVENTINFO_HH
#define BACONANA_DATAFORMATS_TGENEVENTINFO_HH

#include <TObject.h>

namespace baconhep 
{
  class TGenEventInfo : public TObject
  {
    public:
      TGenEventInfo():
      id_1(0),  id_2(0),
      x_1(0),   x_2(0),
      scalePDF(0)
      {}
      ~TGenEventInfo(){}

      int   id_1,  id_2;   // parton flavor PDG ID
      float x_1,   x_2;    // parton momentum fraction
      float scalePDF;      // Q-scale used for PDF evaluation

    ClassDef(TGenEventInfo,1)
  };
}
#endif
