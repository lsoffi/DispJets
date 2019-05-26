#include "TString.h"
#include "analysisDelphes.cpp+"

void AnalysisDelphes(TString indir, TString outdir, TString tcut, TString tres, float lumi)
{
  analysisDelphes ana(indir,outdir,tcut,tres,lumi); 
  ana.run();
}
