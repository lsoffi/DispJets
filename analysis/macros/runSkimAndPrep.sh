#! /bin/sh

path="../"
sample=$1

echo "Running sample: " ${sample}
root -l -b -q "doSkimAndPrep.C("\"${path}\"","\"${sample}\"")"
echo "Done skimming"
