Steps:

1) Make datacards (./runAnalysis.sh) in macro dir.

2) Run combine and extract jsons for limits (combineDispJets.sh) in combo dir (need HiggsCombine and CombineTool setup)
My dir: /afs/cern.ch/work/m/mzientek/private/ggttCombination/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/CombineHarvester/CombineTools/scripts/DispJets/

3) Make 1D plot and rootfile of the limits (python plot1Dlimits.py DISP [MASS] [TRES] [LUMI in pb])
Output root file is stored in same directory as the input jsons.

4) Use the output rootfile and compare with run2 and theorists results (python overlayLimits.py --suffix [ADDNAME] --theory --orig)

