#include "analysisDelphes.hh"

analysisDelphes::analysisDelphes(TString indir, TString outdir, TString t_cut, TString t_res, float setlumi)
{

  // initialize constants
  inpath        = indir;        // input dir 
  outpath       = outdir;       // output dir
  res           = t_res;        // resolution in ps
  tcut          = t_cut;        // jet time cut in ns
  ptcut         = "30";         // jet pt cut in GeV
  lumi          = setlumi*1000; // translate lumi to pb^-1

  // set up cut string
  /*  TString cut0 = Form("(jetPUPPIPT_[0] > %s && jetPUPPITuncorr_[0] > %s)",ptcut.Data(),tcut.Data());
  TString cut1 = Form("(jetPUPPIPT_[1] > %s && jetPUPPITuncorr_[1] > %s)",ptcut.Data(),tcut.Data());
  TString cut2 = Form("(jetPUPPIPT_[2] > %s && jetPUPPITuncorr_[2] > %s)",ptcut.Data(),tcut.Data());
  TString cut3 = Form("(jetPUPPIPT_[3] > %s && jetPUPPITuncorr_[3] > %s)",ptcut.Data(),tcut.Data());
  TString cut4 = Form("(jetPUPPIPT_[4] > %s && jetPUPPITuncorr_[4] > %s)",ptcut.Data(),tcut.Data());
  */
  cut = "nJetsPUPPI0p1_>=2";
  
  std::cout<<cut<<" "<<cut.Data()<<std::endl;
  //  cut  = Form("(%s || %s || %s || %s|| %s)",cut0.Data(),cut1.Data(),cut2.Data(),cut3.Data(),cut4.Data());
  //  cut = Form("Sum$(1*(jetPUPPIPT_>%s && jetPUPPITuncorr_>%s*0.000000001))>=1",ptcut.Data(),tcut.Data());
  // input names
  //s_file.push_back("XXto4Q_M100_CT100mm");
  //s_file.push_back("dispjets");
  //  s_file.push_back("xxqqqq_m50_ct0mm");
  s_file.push_back("ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl10000_ev100000");
  s_file.push_back("ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl1000_ev100000");
  s_file.push_back("ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl500_ev100000");


  //b_file.push_back("QCD");
  nSig = s_file.size();
  nBkg = b_file.size();

  // setup xsec values
  for (int f = 0; f < nSig; f++){
    xsec[s_file[f]] = 0.001;//0.1994; // signal xsec = 0.001pb
  }
  xsec["QCD"] = 10000; 
 
  // output files
  fout = TFile::Open(Form("%sdispjets.root",outpath.Data()),"RECREATE");

}// end analysisDelphes

analysisDelphes::~analysisDelphes()
{

  delete fout;

}// end ~analysisDelphes

void analysisDelphes::run()
{
  double eff[3];
  double ctau[3]={500,1000,10000};
  // get number of events above cut
  std::cout<<nSig<<std::endl;
  for (int f = 0; f < nSig; f++){
    std::cout<< "  SIGNAL : "<<s_file[f]<<std::endl;
    vals[s_file[f]] = applySel(s_file[f]);
    eff[f] = vals[s_file[f]] ;
    vals[s_file[f]] = applyNorm(vals[s_file[f]],xsec[s_file[f]]);
  }

  TGraph* g =new TGraph(3, ctau, eff);
  TCanvas* c = new TCanvas("c","",1);
  c->cd();
  g->Draw("APE");
  c->SaveAs("/afs/cern.ch/user/s/soffi/www/DispJets/LimitsDelphes/eff_vs_ctau.png");


  float sum_bkg_val = 0; 
  for (int f = 0; f < nBkg; f++){
    vals[b_file[f]] = applySel(b_file[f]);
    vals[b_file[f]] = applyNorm(vals[b_file[f]],xsec[b_file[f]]);
    sum_bkg_val += vals[b_file[f]];
  }
  vals["bkg"] = sum_bkg_val;
  vals["bkg"] = 0.00000001; // set bkg to small number for combine

  // write out datacard
  for (int f = 0; f < nSig; f++){
    makeCard(s_file[f]);
  }

  //make plots
 TFile* fplot = TFile::Open(Form("%sdispjets.root",outpath.Data()),"READ");
 TH1F* h_500;
 TH1F* h_1000;
 TH1F* h_10000;
 TLegend* leg = new TLegend(0.6, 0.7, 0.85, 0.9);
 h_500= (TH1F*)fplot->Get("h_njets_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl500_ev100000");
 h_1000= (TH1F*)fplot->Get("h_njets_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl1000_ev100000");
 h_10000= (TH1F*)fplot->Get("h_njets_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl10000_ev100000");

 h_500->SetLineColor(kRed);
 h_1000->SetLineColor(kBlue);
 h_10000->SetLineColor(kGreen);
 leg->AddEntry(h_500, "c#tau= 500mm", "L");
 leg->AddEntry(h_1000, "c#tau= 1000mm", "L");
 leg->AddEntry(h_10000, "c#tau= 10000mm", "L");

 h_500->DrawNormalized();
 leg->Draw("same");
 h_1000->DrawNormalized("same"); 
 h_10000->DrawNormalized("same");
 c->SaveAs("/afs/cern.ch/user/s/soffi/www/DispJets/LimitsDelphes/nJets.png");
 c->SaveAs("/afs/cern.ch/user/s/soffi/www/DispJets/LimitsDelphes/nJets.pdf");

 h_500= (TH1F*)fplot->Get("h_njets0p1_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl500_ev100000");
 h_1000= (TH1F*)fplot->Get("h_njets0p1_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl1000_ev100000");
 h_10000= (TH1F*)fplot->Get("h_njets0p1_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl10000_ev100000");

 h_500->SetLineColor(kRed);
 h_1000->SetLineColor(kBlue);
 h_10000->SetLineColor(kGreen);

 h_500->DrawNormalized();
 leg->Draw("same");
 h_1000->DrawNormalized("same"); 
 h_10000->DrawNormalized("same");
 c->SaveAs("/afs/cern.ch/user/s/soffi/www/DispJets/LimitsDelphes/njets0p1.png");
 c->SaveAs("/afs/cern.ch/user/s/soffi/www/DispJets/LimitsDelphes/njets0p1.pdf");


 h_500= (TH1F*)fplot->Get("h_njets0p5_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl500_ev100000");
 h_1000= (TH1F*)fplot->Get("h_njets0p5_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl1000_ev100000");
 h_10000= (TH1F*)fplot->Get("h_njets0p5_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl10000_ev100000");

 h_500->SetLineColor(kRed);
 h_1000->SetLineColor(kBlue);
 h_10000->SetLineColor(kGreen);

 h_500->DrawNormalized();
 h_1000->DrawNormalized("same"); 
 leg->Draw("same");
 h_10000->DrawNormalized("same");
 c->SaveAs("/afs/cern.ch/user/s/soffi/www/DispJets/LimitsDelphes/njets0p5.png");
 c->SaveAs("/afs/cern.ch/user/s/soffi/www/DispJets/LimitsDelphes/njets0p5.pdf");

 h_500= (TH1F*)fplot->Get("h_njets1p0_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl500_ev100000");
 h_1000= (TH1F*)fplot->Get("h_njets1p0_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl1000_ev100000");
 h_10000= (TH1F*)fplot->Get("h_njets1p0_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl10000_ev100000");

 h_500->SetLineColor(kRed);
 h_1000->SetLineColor(kBlue);
 h_10000->SetLineColor(kGreen);

 h_500->DrawNormalized();
 leg->Draw("same");
 h_1000->DrawNormalized("same"); 
 h_10000->DrawNormalized("same");
 c->SaveAs("/afs/cern.ch/user/s/soffi/www/DispJets/LimitsDelphes/njets1p0.png");
 c->SaveAs("/afs/cern.ch/user/s/soffi/www/DispJets/LimitsDelphes/njets1p0.pdf");

 h_500= (TH1F*)fplot->Get("h_jetT_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl500_ev100000");
 h_1000= (TH1F*)fplot->Get("h_jetT_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl1000_ev100000");
 h_10000= (TH1F*)fplot->Get("h_jetT_ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_mh125_mx50_pl10000_ev100000");

 h_500->SetLineColor(kRed);
 h_1000->SetLineColor(kBlue);
 h_10000->SetLineColor(kGreen);

 h_500->DrawNormalized();
 leg->Draw("same");
 h_1000->DrawNormalized("same"); 
 h_10000->DrawNormalized("same");
 c->SetLogy();
 c->SaveAs("/afs/cern.ch/user/s/soffi/www/DispJets/LimitsDelphes/jetT.png");
 c->SaveAs("/afs/cern.ch/user/s/soffi/www/DispJets/LimitsDelphes/jetT.pdf");



 
}// end run

float analysisDelphes::applySel(TString file)
{
  
  // input file
  TFile *f = TFile::Open(Form("%s%s.root",inpath.Data(),file.Data()),"READ");
  std::cout<<file.Data()<<std::endl;
  if (f == NULL){ std::cout << "File not found: " << inpath << file << std::endl; return -1; }
  // input tree
  TTree *t = (TTree*)f->Get("t1");
  std::cout<<t->GetEntries()<<std::endl;
  if (t == NULL){ std::cout << "Tree not found" << std::endl; return -1; }
  Int_t nJetsPUPPI_;
  Int_t nJetsPUPPI0p1_;
  Int_t nJetsPUPPI0p5_;
  Int_t nJetsPUPPI1p0_;

  std::cout<<t->GetEntries()<<std::endl;

  // get total number (denominator) 
  t->SetBranchAddress("nJetsPUPPI_",&nJetsPUPPI_);
  std::cout<<t->GetEntries()<<std::endl;
  t->SetBranchAddress("nJetsPUPPI0p1_",&nJetsPUPPI0p1_);
  t->SetBranchAddress("nJetsPUPPI0p5_",&nJetsPUPPI0p5_);
  t->SetBranchAddress("nJetsPUPPI1p0_",&nJetsPUPPI1p0_);

  float sum_weight = 0;

  /*  for (unsigned int entry = 0; entry < t->GetEntries(); entry++){
    t->GetEntry(entry);
    sum_weight +=1;
  }
  */  
  std::cout<<"exit from loop" <<std::endl;

  // apply cut
  TH1F * h_num = new TH1F("h_num","",2000,0,100);
  t->Draw("nJetsPUPPI0p1_>> h_num","nJetsPUPPI0p1_>=2");
  TH1F * h_den = (TH1F*)f->Get("nEvents");
  sum_weight=h_den->GetBinContent(1);
  
  TH1F * h_njets = new TH1F(Form("h_njets_%s",file.Data()),"",10,0,10);
  t->Draw(Form("nJetsPUPPI_>> h_njets_%s",file.Data()),"");
  TH1F * h_njets0p1 = new TH1F(Form("h_njets0p1_%s",file.Data()),"",10,0,10);
  t->Draw(Form("nJetsPUPPI0p1_>> h_njets0p1_%s",file.Data()),"");
  TH1F * h_njets0p5 = new TH1F(Form("h_njets0p5_%s",file.Data()),"",10,0,10);
  t->Draw(Form("nJetsPUPPI0p5_>> h_njets0p5_%s",file.Data()),"");
  TH1F * h_njets1p0 = new TH1F(Form("h_njets1p0_%s",file.Data()),"",10,0,10);
  t->Draw(Form("nJetsPUPPI1p0_>> h_njets1p0_%s",file.Data()),"");

  TH1F* h_jetT = new TH1F(Form("h_jetT_%s",file.Data()),"", 100,-2,10);
  t->Draw(Form("jetPUPPITcorr_*100000000>>h_jetT_%s",file.Data()),"");
  // save output histo
  fout->cd();
  //  h->Draw("HIST");
  h_njets->Write();
  h_njets0p1->Write();
  h_njets0p5->Write();
  h_njets1p0->Write();
  h_jetT->Write();
  h_num->Write();
 
  // get number after cuts (numerator) 
  float int_aftercuts = h_num->Integral();
  // calculate efficiency
  float eff = 0;
  if (sum_weight > 0 ) eff = int_aftercuts/sum_weight;
  std::cout<<"return" <<std::endl;
  // delete and finish

  return eff;

}// end applySel

float analysisDelphes::applyNorm(float eff, float xsec)
{
  float numexp = lumi*eff*xsec;
  return numexp;
  std::cout<<numexp<<std::endl;
}// end applyNorm

void analysisDelphes::makeCard(TString sig)
{

  // setup card
  TString cardname = Form("%sdatacard_%s_res%s_lumi%0.0fpb.txt",outpath.Data(),sig.Data(),res.Data(),lumi);
  std::cout << "Writing card " << cardname << std::endl;
  std::ofstream card;
  card.open(cardname); 
  std::cout<<vals[sig]<<std::endl;
  if (card.is_open()){
    card << Form("# Datacard for %s with lumi = %0.1f/pb",sig.Data(),lumi) << endl;
    card << "imax 1 " << endl;
    card << "jmax * " << endl;
    card << "kmax * " << endl; 
    card << "------------------------------------" << endl;
    card << "bin bin1" << endl;
    card << "observation 0" << endl;
    card << "------------------------------------" << endl;
    card << "bin      	bin1		bin1"    << endl;
    card << "process 	sig		bkg"  << endl;
    card << "process	0		1"    << endl;
    card << Form("rate     	%0.3f		%0.13f",vals[sig],vals["bkg"]) << endl;
    card << "------------------------------------" << endl;
    card << "lumi     lnN 	1.03		1.03" << endl;
    card << "trig_eff lnN 	1.50		1.50" << endl;
    card << "sel_eff  lnN 	1.30		1.30" << endl;

  }// end if card is open
  else std::cout << "Unable to open output datacard file" << std::endl;
  // close card
  card.close();

}// end makeCard 
