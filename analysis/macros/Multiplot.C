#include "TString.h"
#include "multiplot.cpp+"

void Multiplot(const TString & path)
{

  std::vector< TString > ctau;
  //ctau.push_back("0mm");
  ctau.push_back("1mm");
  ctau.push_back("10mm");
  ctau.push_back("100mm");
  ctau.push_back("1000mm");
  ctau.push_back("10000mm");

  for (unsigned int i = 0; i < ctau.size(); i++){
    multiplot plot(ctau[i],path);
    plot.go();
  }

}
