# Displaced Jets 

## Set up
```
cmsrel CMSSW_9_2_8
cd CMSSW_9_2_8/src/
cmsenv
git cms-init
git clone git@github.com:mez34/DispJets.git
scram b -j 4
```

## Run analyzer (python dir.)
```
cmsRun configDispJetsAna.py
```
Other files in this directory are for running crab jobs. 
The python configs (crabConfig_) will run crabConfigTemplate.py that reads the data in from the inputfiles_* 

Some notes on the analyzer: 
- Only genjets and genparticles are really needed.
- Loop over genjets is here: https://github.com/mez34/DispJets/blob/master/analysis/plugins/DisplacedJetsAnalyzer.cc#L390-L506
  although this is really not used in plotting / subsequent studies detailed below.  
- More interested in the loop over genparticles: 
  * First get hard-interaction quarks info (ie. these are the basis of the jet) here:
  https://github.com/mez34/DispJets/blob/master/analysis/plugins/DisplacedJetsAnalyzer.cc#L548-L578
  * Store their mom's info here: 
  https://github.com/mez34/DispJets/blob/master/analysis/plugins/DisplacedJetsAnalyzer.cc#L617-L665
  * Loop over genparticles again: https://github.com/mez34/DispJets/blob/master/analysis/plugins/DisplacedJetsAnalyzer.cc#L695
  * Do DeltaR matching to the hard-interaction quarks here: 
  https://github.com/mez34/DispJets/blob/master/analysis/plugins/DisplacedJetsAnalyzer.cc#L729-L767
  * Calculate the distance travelled to the MIP timing layer for each particle:
  https://github.com/mez34/DispJets/blob/master/analysis/plugins/DisplacedJetsAnalyzer.cc#L795-L820

## Make plots (macros dir.)

1) Run skimmer `./runSkimmer.sh [SAMPLE]` or `./runFullSkim.sh` for all samples 
which calls Skimmer.C (which uses skimmer.cpp).

Some notes on skimmer: 
Resolutions are specified here: https://github.com/mez34/DispJets/blob/master/analysis/macros/skimmer.cpp#L360-L365
Skimmer loops over gen particles, depending on which jet they are matched to, adds that particle's time, pT, etc. to the quantities for the jet and appends all this information to the ntuples.

2) Make plots `./runPlotter.sh` which calls Plotter.C (which uses plotter.cpp),
runPlotter.sh specifies directories and Plotter.C has the names of the files. 
This makes the plots for each individual sample.
Unmatched particles (not associated to a jet) are smeared in this script here: https://github.com/mez34/DispJets/blob/master/analysis/macros/plotter.cpp#L151-L155
So if you want to look at their resolutions, it needs to be changed here, not in the skimmer.

3) Overlay plots (from each sample) of one variable using `./runOverlay.sh` which calls Overlay.C (which uses overlay.cpp). For example, this works to compare the jet time as a function of lifetime. 

4) Overlay different variable plots for a single sample using `./runMultiplot.sh` which calls Multiplot.C (which uses multiplot.cpp). For example, this works to compare the time for different time resolutions, but using the same lifetime point. 

5) Run analysis `./runAnalysis.sh` to make datacards for combine. 

## Run combine (combo dir.)

1) Take cards from macro dir (./runAnalysis.sh) as mentioned before.

2) Run combine and extract jsons for limits (combineDispJets.sh) in combo dir (need HiggsCombine and CombineTool setup)

3) Make 1D plot and rootfile of the limits (python plot1Dlimits.py DISP [MASS] [TRES] [LUMI in pb])
Output root file is stored in same directory as the input jsons.

4) Use the output rootfile and compare with run2 and theorists results (python overlayLimits.py --suffix [ADDNAME] --theory --orig)

