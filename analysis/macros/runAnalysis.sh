#! /bin/sh

indir="/eos/cms/store/group/phys_egamma/soffi/XXQQQQ-GEN-SIM/DisplJetAnalyzer/XXQQQQ_m50/"      # input directory
outdir="output_files/"   # output directory
tcut=${1:-"0.1"}         # jet time cut [ns]
tres=${2:-"30"}          # time resolution [0,30,50,70,180,500] allowed
lumi=${3:-"3000"}         # total integrated luminosity [fb]

mkdir -p ${outdir}

echo "Running analysis"
root -l -b -q "Analysis.C("\"${indir}\"","\"${outdir}\"","\"${tcut}\"","\"${tres}\"","${lumi}")"

echo "Done"

echo "Cleaning"
rm *~
rm *_cpp*
rm *_C.d
