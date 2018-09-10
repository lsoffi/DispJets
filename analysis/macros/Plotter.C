#include "TString.h"
#include "plotter.cpp+"

void Plotter(const TString & path, const TString & outpath)
{

  std::vector< TString > ctau;
  ctau.push_back("0mm");
  ctau.push_back("1mm");
  ctau.push_back("10mm");
  ctau.push_back("100mm");
  ctau.push_back("1000mm");
  ctau.push_back("10000mm");

  std::map< TString, TString > infiles;
  std::map< TString, TString > outdirs;
  for (unsigned int i = 0; i < ctau.size(); i++){
    infiles[ctau[i]] = Form("%sntuple_xxqqqq_m50_ct%s_skim.root",path.Data(),ctau[i].Data());
    outdirs[ctau[i]] = Form("%s/XXqqqq_m50_ct%s/",outpath.Data(),ctau[i].Data()); 
    plotter plot(ctau[i],infiles[ctau[i]],outdirs[ctau[i]]);
    plot.go();
  }



}
