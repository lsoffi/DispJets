#! /bin/sh

path="../samples/"
sample=${1:-"ntuple_dispjets"}

echo "Running sample: " ${sample}
root -l -b -q "Skimmer.C("\"${path}\"","\"${sample}\"")"
echo "Done skimming"

echo "Cleaning"
rm *~
rm *_cpp*
rm *_C.d
