#include "analysis.hh"

analysis::analysis(TString indir, TString outdir, TString t_cut, TString t_res, float setlumi)
{

  // initialize constants
  inpath        = indir;        // input dir 
  outpath       = outdir;       // output dir
  res           = t_res;        // resolution in ps
  tcut          = t_cut;        // jet time cut in ns
  ptcut         = "30";         // jet pt cut in GeV
  lumi          = setlumi*1000; // translate lumi to pb^-1

  // set up cut string
  TString cut1 = Form("(jet_pt[0] > %s && jet_smear_%s_t[0] > %s)",ptcut.Data(),res.Data(),tcut.Data());
  TString cut2 = Form("(jet_pt[1] > %s && jet_smear_%s_t[1] > %s)",ptcut.Data(),res.Data(),tcut.Data());
  TString cut3 = Form("(jet_pt[2] > %s && jet_smear_%s_t[2] > %s)",ptcut.Data(),res.Data(),tcut.Data());
  TString cut4 = Form("(jet_pt[3] > %s && jet_smear_%s_t[3] > %s)",ptcut.Data(),res.Data(),tcut.Data());
  cut  = Form("weight*(%s || %s || %s || %s)",cut1.Data(),cut2.Data(),cut3.Data(),cut4.Data());

  // input names
  s_file.push_back("XXto4Q_M100_CT100mm");
  //b_file.push_back("QCD");
  nSig = s_file.size();
  nBkg = b_file.size();

  // setup xsec values
  for (int f = 0; f < nSig; f++){
    xsec[s_file[f]] = 1.0; // signal xsec = 1pb
  }
  xsec["QCD"] = 10000; 
 
  // output files
  fout = TFile::Open(Form("%sdispjets.root",outpath.Data()),"RECREATE");

}// end analysis

analysis::~analysis()
{

  delete fout;

}// end ~analysis

void analysis::run()
{

  // get number of events above cut
  for (int f = 0; f < nSig; f++){
    vals[s_file[f]] = applySel(s_file[f]);
    vals[s_file[f]] = applyNorm(vals[s_file[f]],xsec[s_file[f]]);
  }
  float sum_bkg_val = 0; 
  for (int f = 0; f < nBkg; f++){
    vals[b_file[f]] = applySel(b_file[f]);
    vals[b_file[f]] = applyNorm(vals[b_file[f]],xsec[b_file[f]]);
    sum_bkg_val += vals[b_file[f]];
  }
  vals["bkg"] = sum_bkg_val;

  // write out datacard
  for (int f = 0; f < nSig; f++){
    makeCard(s_file[f]);
  }

}// end run

float analysis::applySel(TString file)
{
  
  // input file
  TFile *f = TFile::Open(Form("%sntuples_%s_skim.root",inpath.Data(),file.Data()),"READ");
  if (f == NULL){ std::cout << "File not found: " << inpath << file << std::endl; return -1; }
  // input tree
  TTree *t = (TTree*)f->Get("dispjets/tree");
  if (t == NULL){ std::cout << "Tree not found" << std::endl; return -1; }
  
  // get total number (denominator) 
  float weight;
  TBranch *b_weight;
  t->SetBranchAddress("weight", &weight, &b_weight);
  float sum_weight = 0;
  for (unsigned int entry = 0; entry < t->GetEntries(); entry++){
    t->GetEntry(entry);
    sum_weight += weight;
  }

  // apply cut
  TH1F * h = new TH1F("h","",150,-5,10);
  if (res=="0") t->Draw("jet_avg_t >> h",Form("%s",cut.Data()));
  else          t->Draw(Form("jet_smear_%s_t >> h",res.Data()),Form("%s",cut.Data()));

  // save output histo
  fout->cd();
  h->Draw("HIST");
  h->Write();
 
  // get number after cuts (numerator) 
  float int_aftercuts = h->Integral();
  // calculate efficiency
  float eff = 0;
  if (sum_weight > 0 ) eff = int_aftercuts/sum_weight;

  // delete and finish
  delete h;
  return eff;

}// end applySel

float analysis::applyNorm(float eff, float xsec)
{
  float numexp = lumi*eff*xsec;
  return numexp;
}// end applyNorm

void analysis::makeCard(TString sig)
{

  // setup card
  TString cardname = Form("%sdatacard_%s_res%s_lumi%0.0fpb.txt",outpath.Data(),sig.Data(),res.Data(),lumi);
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
