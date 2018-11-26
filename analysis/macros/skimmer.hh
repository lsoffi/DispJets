#ifndef __skimmer__
#define __skimmer__

#include "TString.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"
#include <iostream>
#include <vector>

class skimmer
{
public:
  skimmer(TString path, TString sample);
  ~skimmer();
  void  run();
  float calcDeltaT(const float lx, const float bx, const float la, const float ba, const float lo, const float bo);
  float calcAvgT(const float t, const float n); 
  float smearVal(const float res);

private:
  TFile * infile;
  TTree * intree;
  TFile * outfile;
  TTree * outtree;
  TDirectory * outdir;

};

#endif
