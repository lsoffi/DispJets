#include "TString.h"
#include "skimmer.cpp+"

void runSkimmer(const TString & path, const TString & sample)
{
  skimmer skim(path,sample);
  skim.run();

}
