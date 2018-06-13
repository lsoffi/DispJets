#include <iostream>

using namespace std;

float applySel(TString, TString, TString, TFile*, TString, TString);
void makeCard(TString, TString, vector<TString>, map<TString,float>);

void doAnalysis()
{

  TString path = "../";
  TString out  = "output_files/";
  vector<TString> s_file;
  s_file.push_back("XXto4Q_M100_CT100mm");
  int nSig = s_file.size();
  vector<TString> b_file; 
  //b_file.push_back("QCD");
  int nBkg = b_file.size();
  TFile *fout = TFile::Open(Form("%sdispjets.root",out.Data()),"RECREATE");
 
  TString res  = "30";
  TString tcut = "1.0";
  TString ptcut = "30"; 

  TString cut1 = Form("(jet_pt[0] > %s && jet_smear_%s_t[0] > %s)",ptcut.Data(),res.Data(),tcut.Data());
  TString cut2 = Form("(jet_pt[1] > %s && jet_smear_%s_t[1] > %s)",ptcut.Data(),res.Data(),tcut.Data());
  TString cut3 = Form("(jet_pt[2] > %s && jet_smear_%s_t[2] > %s)",ptcut.Data(),res.Data(),tcut.Data());
  TString cut4 = Form("(jet_pt[3] > %s && jet_smear_%s_t[3] > %s)",ptcut.Data(),res.Data(),tcut.Data());
  TString cut  = Form("weight*(%s || %s || %s || %s)",cut1.Data(),cut2.Data(),cut3.Data(),cut4.Data());

  // get number of events above cut
  map<TString, float> vals;
  for (int f = 0; f < nSig; f++){
    vals[s_file[f]] = applySel(path,s_file[f],out,fout,res,cut);
  }
  float sum_bkg_val = 0; 
  for (int f = 0; f < nBkg; f++){
    vals[b_file[f]] = applySel(path,b_file[f],out,fout,res,cut);
    sum_bkg_val += vals[b_file[f]];
  }
  vals["bkg"] = sum_bkg_val;
  // write out datacard
  for (int f = 0; f < nSig; f++){
    makeCard(out,s_file[f],b_file,vals);
  }

}// end doAnalysis

float applySel(TString path, TString file, TString out, TFile* fout, TString res, TString cut){

  // input file
  TFile *f = TFile::Open(Form("%sntuples_%s_skim.root",path.Data(),file.Data()),"READ");
  if (f == NULL){ cout << "File not found: " << path << file << endl; return -1; }
  // input tree
  TTree *t = (TTree*)f->Get("dispjets/tree");
  if (t == NULL){ cout << "Tree not found" << endl; return -1; }

  // apply cut
  TH1F * h = new TH1F("h","",150,-5,10);
  t->Draw(Form("jet_smear_%s_t >> h",res.Data()),Form("%s",cut.Data()));

  fout->cd();
  // save output histo
  h->Draw("HIST");
  h->Write();
  
  float integral = h->Integral();
  return integral;

}// end applySel

void makeCard(TString dir, TString sig, vector<TString> bkgs, map<TString,float> vals){

  // setup card
  TString cardname = Form("%sdatacard_%s.txt",dir.Data(),sig.Data());
  std::cout << "Writing card " << cardname << std::endl;
  std::ofstream card;
  card.open(cardname); 
  
  if (card.is_open()){
    card << Form("# Datacard for %s ",sig.Data()) << endl;
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
