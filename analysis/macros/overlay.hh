#ifndef __overlay__
#define __overlay__

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"

typedef std::map< TString, TString > Samples;

class overlay
{
public:
  overlay(std::vector< TString >, TString, std::map< TString, TString > );
  ~overlay();
  void go();

private:
  std::vector<TString>	ctau;
  Samples		infiles;
  TString		path;

};

#endif
