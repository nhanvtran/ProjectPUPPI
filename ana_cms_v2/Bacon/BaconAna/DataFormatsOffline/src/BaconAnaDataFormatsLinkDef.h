#ifndef BACONANA_DATAFORMATS_LINKDEF_H
#define BACONANA_DATAFORMATS_LINKDEF_H
#include "../interface/TEventInfo.hh"
#include "../interface/TGenEventInfo.hh"
#include "../interface/TGenParticle.hh"
#include "../interface/TElectron.hh"
#include "../interface/TMuon.hh"
#include "../interface/TTau.hh"
#include "../interface/TAddJet.hh"
#include "../interface/TJet.hh"
#include "../interface/TPhoton.hh"
#include "../interface/TVertex.hh"
#include "../interface/TPFPart.hh"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace baconhep;

#pragma link C++ class baconhep::TEventInfo+;
#pragma link C++ class baconhep::TGenEventInfo+;
#pragma link C++ class baconhep::TGenParticle+;
#pragma link C++ class baconhep::TElectron+;
#pragma link C++ class baconhep::TMuon+;
#pragma link C++ class baconhep::TTau+;
#pragma link C++ class baconhep::TJet+;
#pragma link C++ class baconhep::TAddJet+;
#pragma link C++ class baconhep::TPhoton+;
#pragma link C++ class baconhep::TVertex+;
#pragma link C++ class baconhep::TPFPart+;
#endif
