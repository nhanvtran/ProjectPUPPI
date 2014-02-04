# ----------------------------------------
# inputs 

FASTJETLOCATION=/Users/ntran/Research/CMS/PhysicsAnalysis/ProjectPUPPI/V2/fastjet-install/
FASTJETINC :=`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`
FASTJETLIB :=`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins`

ROOTINC := `root-config --glibs --cflags`

CC = c++

# ----------------------------------------


all: puppiAnalyzer

puppiAnalyzer: puppiAnalyzer.cpp
	$(CC) -o puppiAnalyzer $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) -lvectorDict -lJetCleanser puppiAnalyzer.cpp puppiContainer.cc puppiContainer.hh NoTrees.cc NoTrees.hh -I.

clean:
	rm puppiAnalyzer

    