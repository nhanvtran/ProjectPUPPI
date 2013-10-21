void runPuppiRun(int events, int puscenario){
	gROOT->ProcessLine(".x setupFastjet.C");
	gROOT->ProcessLine(".L puppiContainer.cc++");
	gROOT->ProcessLine(".L puppiAnalyzer.cc++");
    //puppiAnalyzer(10,20);
	puppiAnalyzer(events,puscenario);    
}
