#ifndef __overlay__
#define __overlay__

#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"
#include <iostream>

typedef std::map< TString, TFile* > Samples; 
typedef std::vector< TLegend* >     LegVec;
typedef std::vector< TCanvas* >	    CanVec;
typedef std::vector< TString >      Histos;
typedef std::vector< TH1F* >        TH1Vec;
typedef std::vector< TH1Vec >       TH1VecVec;
typedef std::map< TString, TString> Names;
typedef std::map< TString, Color_t> Colors;

class overlay
{
public:
  overlay(std::vector< TString >, TString, Names, Colors);
  ~overlay();
  void go();
  void setupplots();
  void setupextras();
  void prepcanvas();
  void drawplots();
  void saveplots();

private:
  std::vector<TString>	ctau;
  Samples		infile;
  TString		path;
  unsigned int		nfiles;
  TString		odir;
  TFile*		ofile;
  Histos		hnames;
  Histos		hextra; 
  unsigned int		nh;
  TH1VecVec		histos;
  TH1VecVec		ohistos;
  CanVec		canv;
  LegVec		leg;
  Colors		colors;

};

#endif
