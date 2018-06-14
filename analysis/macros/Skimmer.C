#include "TString.h"
#include "skimmer.cpp+"

void Skimmer(const TString & path, const TString & sample)
{
  skimmer skim(path,sample);
  skim.run();

}
