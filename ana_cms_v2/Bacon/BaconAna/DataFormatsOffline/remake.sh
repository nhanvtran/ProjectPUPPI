#!/bin/bash

rm *.o
rm *.so
#for x in ` ls src/*.cc ` ; do g++  -O2 -pipe -Wall -W -Woverloaded-virtual -pthread -m64 -Iinclude -c $x; done
#g++ -dynamiclib -single_module -install_name $PWD/Bacon.so -O2 -mmacosx-version-min=10.8 -m64 `ls *.o` -o  Bacon.so -L $PWD -Wl,-rpath,$PWD  -L $PWD/../src  

rootcint -f src/BaconDict.cc -c `ls interface/*`  src/BaconAnaDataFormatsLinkDef.h
for x in `ls src/*.cc`; do g++   -I. `root-config --cflags` -c $x; done
g++ -shared -single_module -install_name $PWD/libBacon.so -O2 -mmacosx-version-min=10.8 -m64 `root-config --glibs --cflags` `ls *.o`  -o  libBacon.so -L $PWD -Wl,-rpath,$PWD  -L $PWD/src