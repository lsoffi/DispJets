#! /bin/sh

indir="/eos/cms/store/user/pablom/LLPNTuples/ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl10000_ev100000Tuple/"
outdir="output_files_Delphes/"   # output directory
tcut=${1:-""}         # jet time cut [ns]
tres=${2:-"30"}          # time resolution [0,30,50,70,180,500] allowed
lumi=${3:-"3000"}         # total integrated luminosity [fb]

mkdir -p ${outdir}

echo "Running analysisDelphes"
echo "AnalysisDelphes.C("\"${indir}\"","\"${outdir}\"","\"${tcut}\"","\"${tres}\"","${lumi}")"
root -l -b -q "AnalysisDelphes.C("\"${indir}\"","\"${outdir}\"","\"${tcut}\"","\"${tres}\"","${lumi}")"

echo "Done"

echo "Cleaning"
rm *~
rm *_cpp*
rm *_C.d
