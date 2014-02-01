{

    // at cmslpc
//	gSystem->Load("/uscms_data/d2/ntran/physics/Jets/PUPPI/fastjet/myFastjet304/lib/libfastjet.so");
//    gSystem->Load("/uscms_data/d2/ntran/physics/Jets/PUPPI/fastjet/myFastjet304/lib/libfastjettools.so");

    // at home
    gSystem->Load("/Users/ntran/Research/CMS/PhysicsAnalysis/boost2013/fastjet/fastjet-install/lib/libfastjet.so");
    gSystem->Load("/Users/ntran/Research/CMS/PhysicsAnalysis/boost2013/fastjet/fastjet-install/lib/libfastjettools.so");

    gROOT->ProcessLine(".L NoTrees.cc++"); 
    gSystem->Load("NoTrees_cc.so");

    // at cmslpc
//	gSystem->AddIncludePath("-I/uscms_data/d2/ntran/physics/Jets/PUPPI/fastjet/myFastjet304/include");
    // at home
    gSystem->AddIncludePath("/Users/ntran/Research/CMS/PhysicsAnalysis/boost2013/fastjet/fastjet-install/include/");      
}
