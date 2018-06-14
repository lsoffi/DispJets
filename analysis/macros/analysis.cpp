#include "analysis.hh"

analysis::analysis()
{

  // initialize const
  res           = "30";  // resolution in ps
  tcut          = "1.0"; // jet time cut in ns
  ptcut         = "30";  // jet pt cut in GeV
  float setlumi = 2.6;   // lumi in fb^-1
  lumi = setlumi*1000;   // translate lumi to pb^-1

  // set up cut string
  TString cut1 = Form("(jet_pt[0] > %s && jet_smear_%s_t[0] > %s)",ptcut.Data(),res.Data(),tcut.Data());
  TString cut2 = Form("(jet_pt[1] > %s && jet_smear_%s_t[1] > %s)",ptcut.Data(),res.Data(),tcut.Data());
  TString cut3 = Form("(jet_pt[2] > %s && jet_smear_%s_t[2] > %s)",ptcut.Data(),res.Data(),tcut.Data());
  TString cut4 = Form("(jet_pt[3] > %s && jet_smear_%s_t[3] > %s)",ptcut.Data(),res.Data(),tcut.Data());
  cut  = Form("weight*(%s || %s || %s || %s)",cut1.Data(),cut2.Data(),cut3.Data(),cut4.Data());
 
  // output files
  TString path = "../";
  TString out  = "output_files/";
  fout = TFile::Open(Form("%sdispjets.root",out.Data()),"RECREATE");
 

}// end analysis

analysis::~analysis()
{

  delete fout;

}// end ~analysis

void analysis::run()
{

}// end run

float analysis::applySel(TString path, TString file)
{
  
  std::cout << " here " << std::endl;
  float eff = 0;
  return eff;

}// end applySel

float analysis::applyNorm(float eff, float xsec)
{
  float numexp = lumi*eff*xsec;
  //std::cout << "Inputs  : " << eff << " " << xsec << " " << lumi << std::endl;
  //std::cout << "#Events : " << numexp << std::endl; 
  return numexp;
}// end applyNorm

void analysis::makeCard(TString dir, TString sig, std::vector<TString> bkgs, std::map<TString,float> vals)
{

  // setup card
  TString cardname = Form("%sdatacard_%s.txt",dir.Data(),sig.Data());
  std::cout << "Writing card " << cardname << std::endl;
  std::ofstream card;
  card.open(cardname); 
  
  if (card.is_open()){
    card << Form("# Datacard for %s with lumi = %0.1f/pb",sig.Data(),lumi) << endl;
    card << "imax 1 " << endl;
    card << "jmax * " << endl;
    card << "kmax * " << endl; 
    card << "------------------------------------" << endl;
    card << "bin 1" << endl;
    card << "observation 0" << endl;
    card << "------------------------------------" << endl;
    card << "bin      	1		1" << endl;
    card << "process 	sig		bkg" << endl;
    card << Form("rate     	%0.3f		%0.3f",vals[sig],vals["bkg"]) << endl;
    card << "------------------------------------" << endl;
    card << "lumi     	1.03		1.03" << endl;
    card << "trig_eff 	1.50		1.50" << endl;
    card << "sel_eff  	1.30		1.30" << endl;

  }// end if card is open
  else std::cout << "Unable to open output datacard file" << std::endl;
  // close card
  card.close();

}// end makeCard 
