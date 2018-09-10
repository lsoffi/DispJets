#include "TString.h"
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

  std::map< TString, TString > infiles;
  for (unsigned int i = 0; i < ctau.size(); i++){
    infiles[ctau[i]] = Form("%s/XXqqqq_m50_ct%s/plots_dispjets.root",path.Data(),ctau[i].Data()); 
  }

  overlay plot(ctau,path,infiles);
  plot.go();

}
