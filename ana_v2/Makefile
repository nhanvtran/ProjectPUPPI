# ----------------------------------------
# inputs 

FASTJETLOCATION=/Users/ntran/Documents/Research/CMS/PhysicsAnalysis/fastjet/fastjet-install/bin
FASTJETINC :=`$(FASTJETLOCATION)/fastjet-config --cxxflags --plugins`
FASTJETLIB :=`$(FASTJETLOCATION)/fastjet-config --libs --plugins`

ROOTINC := `root-config --glibs --cflags` -lTMVA

CC = g++

# ----------------------------------------


all: puppiAnalyzer puppiJetAnalyzer puppiMetAnalyzer

puppiAnalyzer: puppiAnalyzer.cpp
	$(CC) -o puppiAnalyzer    $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) -lvectorDict -lJetCleanser puppiAnalyzer.cpp    puppiContainer.cc       puppiTMVAContainer.cc  NoTrees.cc  -I.

puppiJetAnalyzer: puppiJetAnalyzer.cpp puppiCleanContainer.cc puppiTMVAContainer.cc
	$(CC) -o puppiJetAnalyzer $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) -lvectorDict -lJetCleanser puppiJetAnalyzer.cpp puppiCleanContainer.cc  puppiTMVAContainer.cc  NoTrees.cc  -I.

puppiMetAnalyzer: puppiMetAnalyzer.cpp puppiCleanContainer.cc puppiTMVAContainer.cc
	$(CC) -o puppiMetAnalyzer $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) -lvectorDict puppiMetAnalyzer.cpp puppiCleanContainer.cc  puppiTMVAContainer.cc  NoTrees.cc  -I.

clean:
	rm puppiAnalyzer
	rm puppiJetAnalyzer
	rm puppiMetAnalyzer    