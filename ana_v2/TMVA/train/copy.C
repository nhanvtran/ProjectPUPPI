void copy(std::string iName="OutputTmp.root",std::string iCut="pu == 0") { 
  TFile *lFile = new TFile(iName.c_str()); 
  TTree *lTree = (TTree*) lFile->FindObjectAny("tree");

  TFile *lFile  = new TFile("Output.root","RECREATE"); 
  TTree *lOTree = (TTree*) lTree->CopyTree(iCut.c_str());
  lOTree->Write();
}
