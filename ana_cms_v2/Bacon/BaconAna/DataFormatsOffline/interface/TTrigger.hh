#ifndef BACONANA_DATAFORMATS_TTRIGGER_HH
#define BACONANA_DATAFORMATS_TTRIGGER_HH

#include "../interface/TriggerRecord.hh"  // class to handle user specified trigger info
#include "../interface/BaconAnaDefs.hh"

namespace baconhep
{
  class TTrigger 
  {
    public:
    TTrigger(std::string iFileName);
      ~TTrigger(){}
    //Methods
    int  getTriggerBit(std::string iName);
    int* getTriggerObjectBits(std::string iName);
    int  getTriggerObjectBit(std::string iName,std::string iObjName);
    bool pass(std::string iName,TriggerBits &iTrig);
    bool passObj(std::string iName,TriggerObjects &iTrigObj);

    std::vector<baconhep::TriggerRecord> fRecords;
  };
}
#endif
