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
  analysis(TString indir, TString outdir, TString t_cut, TString t_res, float setlumi);
  ~analysis();
  void  run();
  float applySel(TString file);
  float applyNorm(float eff, float xsec);
  void  makeCard(TString sig); 

private:
  TString 		res;
  TString 		tcut;
  TString 		ptcut;
  float   		lumi;
  TString 		cut;

  TString 		inpath;
  map<TString,float>	xsec;
  map<TString,float>	vals;
  std::vector<TString>	s_file;
  std::vector<TString>	b_file;
  int			nSig; 
  int			nBkg;
  TString 		outpath;
  TFile * 		fout;
  

};

#endif
