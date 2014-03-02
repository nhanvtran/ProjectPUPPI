#include "../interface/TTrigger.hh"

#include <string.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <limits> 

using namespace baconhep;

TTrigger::TTrigger(std::string iFileName) { 
  std::ifstream lFile(iFileName.c_str());
  assert(lFile.is_open());
  int         lTrigIndex       = -1;
  int         lTrigObjIndex    = -1;
  std::string lLastTrigName = "";
  int         lLastLep      = -1;
  while(lFile.good()) { 
    std::string pTrigName; 
    std::string pFilterName;
    int         pLep;
    lFile >> pTrigName;
    if(pTrigName.find("//") != std::string::npos) {
      lFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }
    lFile >> pFilterName >> pLep;
    if(pTrigName == "") continue;
    if(lLastTrigName != pTrigName) {
      lTrigIndex++;
      fRecords.push_back(baconhep::TriggerRecord(pTrigName,lTrigIndex));
    }
    if(lLastLep      != pLep) lTrigObjIndex++;
    fRecords.back().objectMap.push_back(std::pair<std::string, int>(pFilterName,lTrigObjIndex));
    lLastTrigName = pTrigName;
    lLastLep      = pLep;
    lFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
}
int TTrigger::getTriggerBit(std::string iName) { 
  int lId = -1;
  for(unsigned int i0 = 0; i0 < fRecords.size(); i0++) { 
    if(iName == fRecords[i0].hltPattern) lId = i0;    
  }
  if(lId == -1) std::cout << "=== Missing Trigger ==" << iName << std::endl;
  return lId;
}
int* TTrigger::getTriggerObjectBits(std::string iName) { 
  int lId = getTriggerBit(iName);
  int *lOut;
  if(lId == -1) {lOut = new int[1]; lOut[0] = -1; return lOut;}
  int lNBins = fRecords[lId].objectMap.size()+1;
  lOut = new int[lNBins];
  lOut[0] = lNBins=1;
  for(unsigned int i0 = 0; i0 < fRecords[lId].objectMap.size(); i0++) { 
    lOut[i0+1] = fRecords[lId].objectMap[i0].second;
  }
  return lOut;
}
int TTrigger::getTriggerObjectBit(std::string iName,std::string iObjName) { 
  int lId = getTriggerBit(iName);
  if(lId == -1) return lId;
  for(unsigned int i0 = 0; i0 < fRecords[lId].objectMap.size(); i0++) { 
    if( iObjName != fRecords[lId].objectMap[i0].first) continue;
    lId = i0;
    break;
  }
  return lId;
}
bool TTrigger::pass(std::string iName,TriggerBits &iTrig) {
  int lId = getTriggerBit(iName);
  if(lId == -1) return false;
  return iTrig[lId];
}
bool TTrigger::passObj(std::string iName,TriggerObjects &iTrigObj) {
  int lId = getTriggerBit(iName);
  if(lId == -1) return false;
  bool lPass = false;
  for(unsigned int i0 = 0; i0 < fRecords[lId].objectMap.size(); i0++) { 
    if(iTrigObj[fRecords[lId].objectMap[i0].second]) lPass = true;
  }
  return lPass;
}
