#ifndef BACONANA_NTUPLER_TRIGGERRECORD_HH
#define BACONANA_NTUPLER_TRIGGERRECORD_HH

#include <vector>
#include <string>
#include <utility>

namespace baconhep {

class TriggerRecord
{
public:
  TriggerRecord(const std::string name="", const unsigned int value=0) {
    hltPattern   = name;
    baconTrigBit = value;
    hltPathName  = "";
    hltPathIndex = (unsigned int)-1;
  }
  ~TriggerRecord(){}

  std::string	hltPattern;    // HLT path name/pattern (wildcards allowed: *,?)
  unsigned int  baconTrigBit;  // bacon trigger bit
  std::string	hltPathName;   // HLT path name in trigger menu
  unsigned int  hltPathIndex;  // HLT path index in trigger menu

  // map between trigger object name and bacon trigger object bit
  std::vector<std::pair<std::string,unsigned int> > objectMap;
};

}
#endif
