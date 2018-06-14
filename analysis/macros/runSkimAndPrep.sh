#! /bin/sh

path="../"
sample=${1:-"ntuples_XXto4Q_M100_CT100mm"}

echo "Running sample: " ${sample}
root -l -b -q "runSkimmer.C("\"${path}\"","\"${sample}\"")"
echo "Done skimming"

echo "Cleaning"
rm *~
rm *_cpp*
rm *_C.d
