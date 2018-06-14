#ifndef __analysis__
#define __analysis__

#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

class analysis
{
public:
  analysis();
  ~analysis();
  void  run();
  float applyNorm(float eff, float xsec);
  float applySel(TString path, TString file);
  void  makeCard(TString dir, TString sig, vector<TString> bkgs, map<TString,float> vals); 

private:
  TString res;
  TString tcut;
  TString ptcut;
  float   lumi;
  TString cut;

  TFile * fout;

};

#endif
