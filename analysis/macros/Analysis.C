#include "TString.h"
#include "analysis.cpp+"

void Analysis(TString indir, TString outdir, TString tcut, TString tres, float lumi)
{
  analysis ana(indir,outdir,tcut,tres,lumi); 
  ana.run();
}
