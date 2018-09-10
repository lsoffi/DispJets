#! /bin/sh

indir="../samples/"       # input directory
outdir="output_files/"   # output directory
tcut=${1:-"1.0"}         # jet time cut [ns]
tres=${2:-"30"}          # time resolution [0,30,180] allowed
lumi=${3:-"30000"}       # total integrated luminosity [pb]

mkdir -p ${outdir}

echo "Running analysis"
root -l -b -q "Analysis.C("\"${indir}\"","\"${outdir}\"","\"${tcut}\"","\"${tres}\"","${lumi}")"

echo "Done"

echo "Cleaning"
rm *~
rm *_cpp*
rm *_C.d
