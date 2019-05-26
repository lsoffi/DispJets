#! /bin/sh

#---------------- Setup inputs ---------------------
mass=${1:-"100"}
res=${1:-"30"}
lumi=${2:-"3000000"}

indir="../output_files/"
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
for cTau in 1 10 100 1000 10000;
do

  #card="${indir}datacard_XXto4Q_M${mass}_CT${cTau}mm_res${res}_lumi${lumi}pb.txt"
  card="${indir}datacard_dispjets_ct${cTau}mm_res${res}_lumi${lumi}pb.txt"
  combine -M AsymptoticLimits -m ${cTau} -d ${card}  -n .limit 
#  combineTool.py -M Asymptotic -m ${cTau} -d ${card} --there -n .limit --parallel 18
  #combine -M Asymptotic -m ${mass} -d ${card} 

done

mv higgsCombine.limit.AsymptoticLimits.mH*.root ${outdir1} 

#---------------- Make jsons -------------------------
for cTau in 1 10 100 1000 10000;
do
  name="${outdir2}datacard_dispjets_ct${cTau}mm_res${res}_lumi${lumi}pb.json" #XXto4Q_M${mass}_CT${cTau}_r${res}_L${lumi}.json"
#  combine -M CollectLimits ${outdir1}*.limit.*mH${cTau}.* -o "${name}"
  combineTool.py -M CollectLimits ${outdir1}*.limit.*mH${cTau}.* -o "${name}"
  echo Made json: ${name} 
done
