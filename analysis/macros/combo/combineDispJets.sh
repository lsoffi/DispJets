#! /bin/sh

#---------------- Setup inputs ---------------------
mass=${1:-"100"}
res=${1:-"30"}
lumi=${2:-"2600"}

indir="/afs/cern.ch/work/m/mzientek/private/DispJets/CMSSW_9_2_8/src/DispJets/analysis/macros/output_files/"
outdir1="dispjets_combo/"
outdir2="dispjets_jsons/"
if [ ! -d "$outdir1" ]; then
  mkdir $outdir1
fi
if [ ! -d "$outdir2" ]; then
  mkdir $outdir2
fi

#---------------- Run combine ----------------------

#for cTau in 1 3 10 30 100 300 1000 3000; 
for cTau in 100;
do

  #card="${indir}datacard_XXto4Q_M${mass}_CT${cTau}mm_res${res}_lumi${lumi}pb.txt"
  card="datacard_XXto4Q_M${mass}_CT${cTau}mm_res${res}_lumi${lumi}pb.txt"
  combineTool.py -M Asymptotic -m ${cTau} -d ${card} --there -n .limit --parallel 18

  #combine -M Asymptotic -m ${mass} -d ${card} 

done

mv higgsCombine.limit.Asymptotic.mH*.root ${outdir1} 

#---------------- Make jsons -------------------------
for cTau in 100;
do
  name="${outdir2}XXto4Q_M${mass}_CT${cTau}_r${res}_L${lumi}.json"
  combineTool.py -M CollectLimits ${outdir1}*.limit.*mH${cTau}.* -o "${name}"
  echo Made json: ${name} 
done
