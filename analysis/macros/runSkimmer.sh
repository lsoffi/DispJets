#! /bin/sh

#path="/eos/cms/store/group/phys_egamma/soffi/XXQQQQ-GEN-SIM/DisplJetAnalyzer/XXQQQQ_m50/"
path="../python/"
sample=${1:-"ntuple_dispjets"}

echo "Running sample: " ${sample}
root -l -b -q "Skimmer.C("\"${path}\"","\"${sample}\"")"
echo "Done skimming"

echo "Cleaning"
rm *~
rm *_cpp*
rm *_C.d
