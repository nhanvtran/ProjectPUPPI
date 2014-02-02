{ 
  gSystem->AddIncludePath("-I../../../../ProjectPUPPI-master/fastjet/include");
  gSystem->Load("../../../../ProjectPUPPI-master/fastjet/src/libEvent.so");
  gSystem->Load("../../../../ProjectPUPPI-master/fasjet/tools/libFastJetTools.so");
  //gSystem->Load("NoTrees_cc.so");
  //gSystem->Load("puppiContainer_cc.so");
  gROOT->ProcessLine(".L ../NoTrees.cc+");
  gROOT->ProcessLine(".L ../puppiContainer.cc+");
  gROOT->ProcessLine(".L ../puppiCleanContainer.cc+");
}
