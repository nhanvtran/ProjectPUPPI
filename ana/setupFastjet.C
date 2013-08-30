{

	gSystem->Load("/uscms_data/d2/ntran/physics/Jets/PUPPI/fastjet/myFastjet304/lib/libfastjet.so");
    gSystem->Load("/uscms_data/d2/ntran/physics/Jets/PUPPI/fastjet/myFastjet304/lib/libfastjettools.so");
    gSystem->Load("NoTrees_cc.so");
    
	gSystem->AddIncludePath("-I/uscms_data/d2/ntran/physics/Jets/PUPPI/fastjet/myFastjet304/include");
	
}
