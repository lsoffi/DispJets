#ifndef __multiplot__
#define __multiplot__

#include "TFile.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TColor.h"
#include <iostream>

typedef std::vector< TString >      TStrVec;
typedef std::vector< TStrVec >      TStrVecVec;
typedef std::vector< TH1F* >        TH1FVec;
typedef std::vector< TH1FVec >      TH1FVecVec;
typedef std::vector< TLegend* >     LegVec;
typedef std::vector< TCanvas* >	    CanVec;
typedef std::map< TString, Color_t> ColMap;

class multiplot
{
public:
  multiplot(TString,TString);
  ~multiplot();
  void go(); 
  void setuphistos();
  void prepcanvas();
  void drawplots();
  void saveplots();
 
private:
  TString		ctau;
  TString		path;
  TString		odir;
  TFile*		fin;
  TStrVecVec		hist;
  TH1FVecVec		h;
  unsigned int		nh;
  std::vector<UInt_t>	nsub;
  TFile*        	ofile;
  LegVec		leg;
  CanVec		canv; 
  ColMap		colors;
 
};
#endif
