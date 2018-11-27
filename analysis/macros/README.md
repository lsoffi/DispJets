# Basic instructions

1) Run skimmer `./runSkimmer.sh [SAMPLE]` or `./runFullSkim.sh` for all samples 
which calls Skimmer.C (which uses skimmer.cpp).

2) Make plots `./runPlotter.sh` which calls Plotter.C (which uses plotter.cpp),
runPlotter.sh specifies directories and Plotter.C has the names of the files. 

3) Overlay plots (from each sample) of one variable using `./runOverlay.sh`

4) Overlay different variable plots for a single sample using `./runMultiplot.sh`

5) Run analysis `./runAnalysis.sh` to make datacards for combine 

See next README for limit plotting instructions
(https://github.com/mez34/DispJets/blob/master/analysis/combo/README.md)
