#include "TString.h"
#include "TColor.h"
#include <iostream>
#include "overlay.hh"
#include "overlay.cpp+"

void Overlay(const TString & path)
{

  std::vector< TString > ctau;
  //ctau.push_back("0mm");
  ctau.push_back("1mm");
  ctau.push_back("10mm");
  ctau.push_back("100mm");
  ctau.push_back("1000mm");
  ctau.push_back("10000mm");

  std::map<TString,Color_t> colors;
  //  colors["0mm"] 	= kBlack;
  colors["1mm"] 	= kMagenta;
  colors["10mm"]	= kCyan+1;
  colors["100mm"]	= kGreen-6;
  colors["1000mm"]	= kAzure+7;
  colors["10000mm"]	= kViolet+2; 

  //std::map< TString, TString > infiles;
  Names infiles;
  for (unsigned int i = 0; i < ctau.size(); i++){
    infiles[ctau[i]] = Form("%s/XXqqqq_m50_ct%s/plots_dispjets.root",path.Data(),ctau[i].Data());
  }

  overlay plot(ctau,path,infiles,colors);
  plot.go();

}
