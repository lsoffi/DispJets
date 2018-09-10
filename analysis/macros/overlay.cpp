#include "overlay.hh"

overlay::overlay(std::vector< TString > cTau, TString inpath, std::map< TString, TString > files)
{

  ctau = cTau;
  path = inpath;
  infiles = files; 

} // end overlay constructor

overlay::~overlay()
{

}// end overlay destructor

void overlay::go()
{

}// end go
