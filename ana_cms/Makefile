# ----------------------------------------
# inputs 

FASTJETLOCATION=/Users/Phil/Physics/JWJ/git_v2/ProjectPUPPI/fastjet-3.0.6
FASTJETINC :=`$(FASTJETLOCATION)/fastjet-config --cxxflags --plugins`
FASTJETLIB :=`$(FASTJETLOCATION)/fastjet-config --libs --plugins`
BACONLIB   :=-L/Users/Phil/Physics/JWJ/git_v2/ProjectPUPPI/ana_cms/Bacon/BaconAna/DataFormatsOffline/ -lBacon  -Wl,-rpath,/Users/Phil/Physics/JWJ/git_v2/ProjectPUPPI/ana_cms/Bacon/BaconAna/DataFormatsOffline/ -lm -ldl
ROOTINC := `root-config --glibs --cflags` -lTMVA

CC = g++

# ----------------------------------------


all: runBacon puppiCMSAnalyzer puppiJetAnalyzer 

puppiJetAnalyzer: puppiJetAnalyzer.cpp  puppiTMVAContainer.cc
	$(CC) -o puppiJetAnalyzer $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) -lvectorDict -lJetCleanser puppiJetAnalyzer.cpp RecoObj.hh   puppiTMVAContainer.cc puppiTMVAContainer.hh puppiCleanContainer.cc puppiCleanContainer.hh NoTrees.cc NoTrees.hh -I.

runBacon: runBacon.cpp
	$(CC) -o runBacon    $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) $(BACONLIB) -lvectorDict -lJetCleanser runBacon.cpp    GenMatcher.cc  -I.

puppiCMSAnalyzer: 
	$(CC) -o puppiCMSAnalyzer $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) $(BACONLIB) -lvectorDict -lJetCleanser puppiCMSAnalyzer.cpp RecoObj.hh puppiContainer.cc puppiContainer.hh NoTrees.cc NoTrees.hh  -I.


clean:
	rm runBacon
	rm puppiCMSAnalyzer
	rm puppiJetAnalyzer