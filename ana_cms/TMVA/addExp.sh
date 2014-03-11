#!/bin/bash

file=$1
root -b -q apply.C+\(\"$file\"\)
mv Output.root BDTOutput.root
root -b -q applypvalue.C+
mv Output.root $file