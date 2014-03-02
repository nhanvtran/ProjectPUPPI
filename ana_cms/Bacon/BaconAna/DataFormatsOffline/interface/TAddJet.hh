#ifndef BACONANA_DATAFORMATS_TADDJET_HH
#define BACONANA_DATAFORMATS_TADDJET_HH

#include <TObject.h>

namespace baconhep
{
  class TAddJet : public TObject
  {
    public:
      TAddJet():
	index(-1),
	pt_p1(0),ptraw_p1(0),eta_p1(0),phi_p1(0),mass_p1(0),area_p1(0),
        pt_p2(0),ptraw_p2(0),eta_p2(0),phi_p2(0),mass_p2(0),area_p2(0),
        pt_t1(0),ptraw_t1(0),eta_t1(0),phi_t1(0),mass_t1(0),area_t1(0),
        pt_t2(0),ptraw_t2(0),eta_t2(0),phi_t2(0),mass_t2(0),area_t2(0),
        pt_t3(0),ptraw_t3(0),eta_t3(0),phi_t3(0),mass_t3(0),area_t3(0),
        pt_t4(0),ptraw_t4(0),eta_t4(0),phi_t4(0),mass_t4(0),area_t4(0),
        pt_f1(0),ptraw_f1(0),eta_f1(0),phi_f1(0),mass_f1(0),area_f1(0),
        pt_f2(0),ptraw_f2(0),eta_f2(0),phi_f2(0),mass_f2(0),area_f2(0),
        c2_0(0),c2_0P2(0),c2_0P5(0),c2_1P0(0),c2_2P0(0)
    {}
      ~TAddJet(){}
    
    unsigned int index;                                    // Index of original jet collection
    float pt_p1,ptraw_p1,eta_p1,phi_p1,mass_p1,area_p1;    // Jet Info with pruning default   
    float pt_p2,ptraw_p2,eta_p2,phi_p2,mass_p2,area_p2;    // Jet Info with pruning Tight 
    float pt_t1,ptraw_t1,eta_t1,phi_t1,mass_t1,area_t1;    // Jet Info with trimming default   
    float pt_t2,ptraw_t2,eta_t2,phi_t2,mass_t2,area_t2;    // Jet Info with trimming Tight 
    float pt_t3,ptraw_t3,eta_t3,phi_t3,mass_t3,area_t3;    // Jet Info with trimming Tightish
    float pt_t4,ptraw_t4,eta_t4,phi_t4,mass_t4,area_t4;    // Jet Info with trimming Small cone
    float pt_f1,ptraw_f1,eta_f1,phi_f1,mass_f1,area_f1;    // Jet Info with filtering default   
    float pt_f2,ptraw_f2,eta_f2,phi_f2,mass_f2,area_f2;    // Jet Info with fitlering Tight 
    float c2_0,c2_0P2,c2_0P5,c2_1P0,c2_2P0;                // Correlation function of varying exponent
    ClassDef(TAddJet,1)
  };
}
#endif
