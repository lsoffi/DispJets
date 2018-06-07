#include <iostream>

// declare types
typedef std::map<TString, TH1F*> TH1map;
typedef std::map<TString, TH2F*> TH2map;

// declare other functions
void run( TString, TString, TFile* );
void histos( TH1map & , TH2map & );
void save1Dplots( TString, TFile* , const TH1map & );
void save2Dplots( TString, TFile* , const TH2map & );
TH1F * MakeTH1FPlot( const TString, const TString, const int, const double, const double, const TString, const TString );
TH2F * MakeTH2FPlot( const TString, const TString, const int, const double, const double, const int, const double, const double, const TString, const TString );
float calcDeltaT(const float, const float, const float, const float, const float, const float);
float calcAvgT(const float, const int);

// ---------------------------------- 
// ----------- START ----------------
// ---------------------------------- 
void getTimingOfJet()
{

  TString path = "../";
  TString out  = "~/www/Plots/DispJets/GenLevelPlots/Timing/";
  std::vector< TString > file;
  file.push_back(Form("%sntuple_dispjets.root",path.Data()));
  int nsamples = file.size();

  TFile *fout = TFile::Open(Form("%splots_dispjets.root",out.Data()),"RECREATE");
 
  for (int f = 0; f < nsamples; f++){
    run(file[f],out,fout);
  } 

}

void histos( TH1map & map , TH2map & map2){

  map["LL_beta"]	= MakeTH1FPlot("LL_beta","",100,0,1,"LL particle #beta","");
  map["LL_cTau"]	= MakeTH1FPlot("LL_cTau","",50,0,50,"LL particle c#tau [cm]","");
  map["nconst"]		= MakeTH1FPlot("nconst","",50,0,50,"Num. jet constituents","");
  map["const_t_ex"]	= MakeTH1FPlot("const_t_ex","",150,-5,10,"Jet constituent time [ns] Example","");
  map["const_t"]	= MakeTH1FPlot("const_t","",150,-5,10,"Jet constituent time [ns]",""); 
  map["const_id"]	= MakeTH1FPlot("const_id","",400,-200,200,"Jet constituent pdgId","");
  map["const_pt"]	= MakeTH1FPlot("const_pt","",100,0,100,"Jet constituent p_T [GeV]","");
  map["const_eta"]	= MakeTH1FPlot("const_eta","",40,0,2.0,"Jet constituent |#eta|","");
  map["jet_t"]		= MakeTH1FPlot("jet_t","",150,-5,10,"Jet time [ns]",""); 
  map["unmatch_t"]	= MakeTH1FPlot("unmatch_t","",150,-5,10,"Unmatched particle time [ns]","");
  map["max_jet_t"]	= MakeTH1FPlot("max_jet_t","",150,-5,10,"Max jet time [ns]","");

  map2["tX1_tX2"]	= MakeTH2FPlot("tX1_tX2","",150,-5,10,150,-5,10,"Time jets from X1 [ns]","Time jets from X2 [ns]");
  map2["tnear_tfar"]	= MakeTH2FPlot("tnear_tfar","",150,-5,10,250,-5,20,"Time close jets [ns]","Time far jets [ns]");

}


TH1F * MakeTH1FPlot(const TString hname, const TString htitle, const int nbins, const double xlow, const double xhigh, const TString xtitle, const TString ytitle){
  TString ytitleNew;
  Float_t binwidth = (xhigh-xlow)/nbins;
  if (ytitle=="") ytitleNew = "Events"; //Form("Events / %2.1f GeV",binwidth);
  else ytitleNew = ytitle;
 
  TH1F * hist = new TH1F(hname.Data(),htitle.Data(),nbins,xlow,xhigh);
  hist->GetXaxis()->SetTitle(xtitle.Data());
  hist->GetYaxis()->SetTitle(ytitleNew.Data());
  hist->Sumw2();
  gStyle->SetOptStat(1111111);
  return hist;
}// end MakeTH1FPlot

TH2F * MakeTH2FPlot(const TString hname, const TString htitle, const Int_t xnbins, const Double_t xlow, const Double_t xhigh, const Int_t ynbins, const Double_t ylow, const Double_t yhigh, const TString xtitle, const TString ytitle){
  TH2F * hist = new TH2F(hname.Data(),htitle.Data(),xnbins,xlow,xhigh,ynbins,ylow,yhigh);
  hist->GetXaxis()->SetTitle(xtitle.Data());
  hist->GetYaxis()->SetTitle(ytitle.Data());
  return hist;
}// end MakeTH2FPlot

void save1Dplots(TString odir, TFile* fout, const TH1map & map){
  fout->cd();
  std::cout << "Saving" << std::endl; 

  for (const auto & h:map){ // loop over histos
    const auto & name = h.first;
    const auto & hist = h.second;

    // save to output file
    hist->Write(hist->GetName(),TObject::kWriteDelete);

    // draw to canvas
    TCanvas *c = new TCanvas("c","c");
    c->cd();
    hist->Draw("HIST");
    c->SetLogy(0);
    c->SaveAs(Form("%s%s.png",odir.Data(),name.Data()));
    c->SaveAs(Form("%s%s.pdf",odir.Data(),name.Data()));
    c->SetLogy(1);
    c->SaveAs(Form("%s%s_log.png",odir.Data(),name.Data()));
    c->SaveAs(Form("%s%s_log.pdf",odir.Data(),name.Data()));

    delete c;
  }// end loop over histos

}// end save1Dplots

void save2Dplots(TString odir, TFile* fout, const TH2map & map){
  fout->cd();
  std::cout << "Saving" << std::endl; 

  for (const auto & h:map){ // loop over histos
    const auto & name = h.first;
    const auto & hist = h.second;

    // save to output file
    hist->Write(hist->GetName(),TObject::kWriteDelete);

    // draw to canvas
    TCanvas *c = new TCanvas("c","c");
    c->cd();
    hist->Draw("HIST");
    c->SaveAs(Form("%s%s.png",odir.Data(),name.Data()));
    c->SaveAs(Form("%s%s.pdf",odir.Data(),name.Data()));

    delete c;
  }// end loop over histos

}// end save2Dplots


void run(TString file, TString out, TFile* fout){

  // get input file
  TFile *f = TFile::Open(Form("%s",file.Data()),"READ");
  if (f == NULL) std::cout << "File not found: " << file << std::endl;

  // get input tree
  TTree *t = (TTree*)f->Get("dispjets/tree");
  if (t == NULL) std::cout << "Tree not found" << std::endl;

  // declaration of leaf types

   Int_t           sample;
   Int_t           run;
   Int_t           lumi;
   Long64_t        event;
   Int_t           ngenjets;
   vector<float>   *genjet_pt;
   vector<float>   *genjet_e;
   vector<float>   *genjet_eta;
   vector<float>   *genjet_phi;
   vector<int>     *genjet_ndaug;
   vector<int>     *genjet_nconst;
   vector<float>   *genjet_vx;
   vector<float>   *genjet_vy;
   vector<float>   *genjet_vz;
   vector<float>   *genjet_i;
   vector<int>     *genjet_match;
   vector<float>   *genjet0_const_st;
   vector<float>   *genjet0_const_id;
   vector<float>   *genjet0_const_pt;
   vector<float>   *genjet0_const_vx;
   vector<float>   *genjet0_const_vy;
   vector<float>   *genjet0_const_vz;
   vector<float>   *genjet1_const_st;
   vector<float>   *genjet1_const_id;
   vector<float>   *genjet1_const_pt;
   vector<float>   *genjet1_const_vx;
   vector<float>   *genjet1_const_vy;
   vector<float>   *genjet1_const_vz;
   vector<float>   *genjet2_const_st;
   vector<float>   *genjet2_const_id;
   vector<float>   *genjet2_const_pt;
   vector<float>   *genjet2_const_vx;
   vector<float>   *genjet2_const_vy;
   vector<float>   *genjet2_const_vz;
   vector<float>   *genjet3_const_st;
   vector<float>   *genjet3_const_id;
   vector<float>   *genjet3_const_pt;
   vector<float>   *genjet3_const_vx;
   vector<float>   *genjet3_const_vy;
   vector<float>   *genjet3_const_vz;
   Int_t           ngenpart;
   vector<int>     *genpar_id;
   vector<int>     *genpar_stat;
   vector<float>   *genpar_pt;
   vector<float>   *genpar_eta;
   vector<float>   *genpar_phi;
   vector<float>   *genpar_vx;
   vector<float>   *genpar_vy;
   vector<float>   *genpar_vz;
   vector<float>   *genpar_Lxy;
   vector<float>   *genpar_Lxyz;
   vector<int>     *genpar_qnum;
   vector<int>     *genpar_match_q0;
   vector<int>     *genpar_match_q1;
   vector<int>     *genpar_match_q2;
   vector<int>     *genpar_match_q3;
   vector<int>     *genpar_match_q4;
   vector<float>   *genpar_lo;
   vector<float>   *genpar_la;
   Int_t           *nmothers;
   vector<int>     *mom_id;
   vector<int>     *mom_stat;
   vector<float>   *mom_e;
   vector<float>   *mom_m;
   vector<float>   *mom_p;
   vector<float>   *mom_pt;
   vector<float>   *mom_eta;
   vector<float>   *mom_phi;
   vector<float>   *mom_vx;
   vector<float>   *mom_vy;
   vector<float>   *mom_vz;
   vector<float>   *mom_beta;
   vector<float>   *mom_gama;
   vector<float>   *mom_Lxy;
   vector<float>   *mom_Lz;
   vector<float>   *mom_Lxyz;
   vector<float>   *mom_ctau;
   vector<float>   *mom_dupl;

   genjet_pt = 0;
   genjet_e = 0;
   genjet_eta = 0;
   genjet_phi = 0;
   genjet_ndaug = 0;
   genjet_nconst = 0;
   genjet_vx = 0;
   genjet_vy = 0;
   genjet_vz = 0;
   genjet_i = 0;
   genjet_match = 0;
   genjet0_const_st = 0;
   genjet0_const_id = 0;
   genjet0_const_pt = 0;
   genjet0_const_vx = 0;
   genjet0_const_vy = 0;
   genjet0_const_vz = 0;
   genjet1_const_st = 0;
   genjet1_const_id = 0;
   genjet1_const_pt = 0;
   genjet1_const_vx = 0;
   genjet1_const_vy = 0;
   genjet1_const_vz = 0;
   genjet2_const_st = 0;
   genjet2_const_id = 0;
   genjet2_const_pt = 0;
   genjet2_const_vx = 0;
   genjet2_const_vy = 0;
   genjet2_const_vz = 0;
   genjet3_const_st = 0;
   genjet3_const_id = 0;
   genjet3_const_pt = 0;
   genjet3_const_vx = 0;
   genjet3_const_vy = 0;
   genjet3_const_vz = 0;
   genpar_id = 0;
   genpar_stat = 0;
   genpar_pt = 0;
   genpar_eta = 0;
   genpar_phi = 0;
   genpar_vx = 0;
   genpar_vy = 0;
   genpar_vz = 0;
   genpar_Lxy = 0;
   genpar_Lxyz = 0;
   genpar_qnum = 0;
   genpar_match_q0 = 0;
   genpar_match_q1 = 0;
   genpar_match_q2 = 0;
   genpar_match_q3 = 0;
   genpar_match_q4 = 0;
   genpar_lo = 0;
   genpar_la = 0;
   mom_id = 0;
   mom_stat = 0;
   mom_e = 0;
   mom_m = 0;
   mom_p = 0;
   mom_pt = 0;
   mom_eta = 0;
   mom_phi = 0;
   mom_vx = 0;
   mom_vy = 0;
   mom_vz = 0;
   mom_beta = 0;
   mom_gama = 0;
   mom_Lxy = 0;
   mom_Lz = 0;
   mom_Lxyz = 0;
   mom_ctau = 0;
   mom_dupl = 0;

   // List of branches
   TBranch        *b_sample;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_ngenjets;   //!
   TBranch        *b_genjet_pt;   //!
   TBranch        *b_genjet_e;   //!
   TBranch        *b_genjet_eta;   //!
   TBranch        *b_genjet_phi;   //!
   TBranch        *b_genjet_ndaug;   //!
   TBranch        *b_genjet_nconst;   //!
   TBranch        *b_genjet_vx;   //!
   TBranch        *b_genjet_vy;   //!
   TBranch        *b_genjet_vz;   //!
   TBranch        *b_genjet_i;   //!
   TBranch        *b_genjet_match;   //!
   TBranch        *b_genjet0_const_st;   //!
   TBranch        *b_genjet0_const_id;   //!
   TBranch        *b_genjet0_const_pt;   //!
   TBranch        *b_genjet0_const_vx;   //!
   TBranch        *b_genjet0_const_vy;   //!
   TBranch        *b_genjet0_const_vz;   //!
   TBranch        *b_genjet1_const_st;   //!
   TBranch        *b_genjet1_const_id;   //!
   TBranch        *b_genjet1_const_pt;   //!
   TBranch        *b_genjet1_const_vx;   //!
   TBranch        *b_genjet1_const_vy;   //!
   TBranch        *b_genjet1_const_vz;   //!
   TBranch        *b_genjet2_const_st;   //!
   TBranch        *b_genjet2_const_id;   //!
   TBranch        *b_genjet2_const_pt;   //!
   TBranch        *b_genjet2_const_vx;   //!
   TBranch        *b_genjet2_const_vy;   //!
   TBranch        *b_genjet2_const_vz;   //!
   TBranch        *b_genjet3_const_st;   //!
   TBranch        *b_genjet3_const_id;   //!
   TBranch        *b_genjet3_const_pt;   //!
   TBranch        *b_genjet3_const_vx;   //!
   TBranch        *b_genjet3_const_vy;   //!
   TBranch        *b_genjet3_const_vz;   //!
   TBranch        *b_ngenpart;   //!
   TBranch        *b_genpar_id;   //!
   TBranch        *b_genpar_stat;   //!
   TBranch        *b_genpar_pt;   //!
   TBranch        *b_genpar_eta;   //!
   TBranch        *b_genpar_phi;   //!
   TBranch        *b_genpar_vx;   //!
   TBranch        *b_genpar_vy;   //!
   TBranch        *b_genpar_vz;   //!
   TBranch        *b_genpar_Lxy;   //!
   TBranch        *b_genpar_Lxyz;   //!
   TBranch        *b_genpar_qnum;   //!
   TBranch        *b_genpar_match_q0;   //!
   TBranch        *b_genpar_match_q1;   //!
   TBranch        *b_genpar_match_q2;   //!
   TBranch        *b_genpar_match_q3;   //!
   TBranch        *b_genpar_match_q4;   //!
   TBranch        *b_genpar_lo;   //!
   TBranch        *b_genpar_la;   //!
   TBranch        *b_nmothers;   //!
   TBranch        *b_mom_id;   //!
   TBranch        *b_mom_stat;   //!
   TBranch        *b_mom_e;   //!
   TBranch        *b_mom_m;   //!
   TBranch        *b_mom_p;   //!
   TBranch        *b_mom_pt;   //!
   TBranch        *b_mom_eta;   //!
   TBranch        *b_mom_phi;   //!
   TBranch        *b_mom_vx;   //!
   TBranch        *b_mom_vy;   //!
   TBranch        *b_mom_vz;   //!
   TBranch        *b_mom_beta;   //!
   TBranch        *b_mom_gama;   //!
   TBranch        *b_mom_Lxy;   //!
   TBranch        *b_mom_Lz;   //!
   TBranch        *b_mom_Lxyz;   //!
   TBranch        *b_mom_ctau;   //!
   TBranch        *b_mom_dupl;   //!

   // set branches
   t->SetBranchAddress("sample", &sample, &b_sample);
   t->SetBranchAddress("run", &run, &b_run);
   t->SetBranchAddress("lumi", &lumi, &b_lumi);
   t->SetBranchAddress("event", &event, &b_event);
   t->SetBranchAddress("ngenjets", &ngenjets, &b_ngenjets);
   t->SetBranchAddress("genjet_pt", &genjet_pt, &b_genjet_pt);
   t->SetBranchAddress("genjet_e", &genjet_e, &b_genjet_e);
   t->SetBranchAddress("genjet_eta", &genjet_eta, &b_genjet_eta);
   t->SetBranchAddress("genjet_phi", &genjet_phi, &b_genjet_phi);
   t->SetBranchAddress("genjet_ndaug", &genjet_ndaug, &b_genjet_ndaug);
   t->SetBranchAddress("genjet_nconst", &genjet_nconst, &b_genjet_nconst);
   t->SetBranchAddress("genjet_vx", &genjet_vx, &b_genjet_vx);
   t->SetBranchAddress("genjet_vy", &genjet_vy, &b_genjet_vy);
   t->SetBranchAddress("genjet_vz", &genjet_vz, &b_genjet_vz);
   t->SetBranchAddress("genjet_i", &genjet_i, &b_genjet_i);
   t->SetBranchAddress("genjet_match", &genjet_match, &b_genjet_match);
   t->SetBranchAddress("genjet0_const_st", &genjet0_const_st, &b_genjet0_const_st);
   t->SetBranchAddress("genjet0_const_id", &genjet0_const_id, &b_genjet0_const_id);
   t->SetBranchAddress("genjet0_const_pt", &genjet0_const_pt, &b_genjet0_const_pt);
   t->SetBranchAddress("genjet0_const_vx", &genjet0_const_vx, &b_genjet0_const_vx);
   t->SetBranchAddress("genjet0_const_vy", &genjet0_const_vy, &b_genjet0_const_vy);
   t->SetBranchAddress("genjet0_const_vz", &genjet0_const_vz, &b_genjet0_const_vz);
   t->SetBranchAddress("genjet1_const_st", &genjet1_const_st, &b_genjet1_const_st);
   t->SetBranchAddress("genjet1_const_id", &genjet1_const_id, &b_genjet1_const_id);
   t->SetBranchAddress("genjet1_const_pt", &genjet1_const_pt, &b_genjet1_const_pt);
   t->SetBranchAddress("genjet1_const_vx", &genjet1_const_vx, &b_genjet1_const_vx);
   t->SetBranchAddress("genjet1_const_vy", &genjet1_const_vy, &b_genjet1_const_vy);
   t->SetBranchAddress("genjet1_const_vz", &genjet1_const_vz, &b_genjet1_const_vz);
   t->SetBranchAddress("genjet2_const_st", &genjet2_const_st, &b_genjet2_const_st);
   t->SetBranchAddress("genjet2_const_id", &genjet2_const_id, &b_genjet2_const_id);
   t->SetBranchAddress("genjet2_const_pt", &genjet2_const_pt, &b_genjet2_const_pt);
   t->SetBranchAddress("genjet2_const_vx", &genjet2_const_vx, &b_genjet2_const_vx);
   t->SetBranchAddress("genjet2_const_vy", &genjet2_const_vy, &b_genjet2_const_vy);
   t->SetBranchAddress("genjet2_const_vz", &genjet2_const_vz, &b_genjet2_const_vz);
   t->SetBranchAddress("genjet3_const_st", &genjet3_const_st, &b_genjet3_const_st);
   t->SetBranchAddress("genjet3_const_id", &genjet3_const_id, &b_genjet3_const_id);
   t->SetBranchAddress("genjet3_const_pt", &genjet3_const_pt, &b_genjet3_const_pt);
   t->SetBranchAddress("genjet3_const_vx", &genjet3_const_vx, &b_genjet3_const_vx);
   t->SetBranchAddress("genjet3_const_vy", &genjet3_const_vy, &b_genjet3_const_vy);
   t->SetBranchAddress("genjet3_const_vz", &genjet3_const_vz, &b_genjet3_const_vz);
   t->SetBranchAddress("ngenpart", &ngenpart, &b_ngenpart);
   t->SetBranchAddress("genpar_id", &genpar_id, &b_genpar_id);
   t->SetBranchAddress("genpar_stat", &genpar_stat, &b_genpar_stat);
   t->SetBranchAddress("genpar_pt", &genpar_pt, &b_genpar_pt);
   t->SetBranchAddress("genpar_eta", &genpar_eta, &b_genpar_eta);
   t->SetBranchAddress("genpar_phi", &genpar_phi, &b_genpar_phi);
   t->SetBranchAddress("genpar_vx", &genpar_vx, &b_genpar_vx);
   t->SetBranchAddress("genpar_vy", &genpar_vy, &b_genpar_vy);
   t->SetBranchAddress("genpar_vz", &genpar_vz, &b_genpar_vz);
   t->SetBranchAddress("genpar_Lxy", &genpar_Lxy, &b_genpar_Lxy);
   t->SetBranchAddress("genpar_Lxyz", &genpar_Lxyz, &b_genpar_Lxyz);
   t->SetBranchAddress("genpar_qnum", &genpar_qnum, &b_genpar_qnum);
   t->SetBranchAddress("genpar_match_q0", &genpar_match_q0, &b_genpar_match_q0);
   t->SetBranchAddress("genpar_match_q1", &genpar_match_q1, &b_genpar_match_q1);
   t->SetBranchAddress("genpar_match_q2", &genpar_match_q2, &b_genpar_match_q2);
   t->SetBranchAddress("genpar_match_q3", &genpar_match_q3, &b_genpar_match_q3);
   t->SetBranchAddress("genpar_match_q4", &genpar_match_q4, &b_genpar_match_q4);
   t->SetBranchAddress("genpar_lo", &genpar_lo, &b_genpar_lo);
   t->SetBranchAddress("genpar_la", &genpar_la, &b_genpar_la);
   t->SetBranchAddress("nmothers", &nmothers, &b_nmothers);
   t->SetBranchAddress("mom_id", &mom_id, &b_mom_id);
   t->SetBranchAddress("mom_stat", &mom_stat, &b_mom_stat);
   t->SetBranchAddress("mom_e", &mom_e, &b_mom_e);
   t->SetBranchAddress("mom_m", &mom_m, &b_mom_m);
   t->SetBranchAddress("mom_p", &mom_p, &b_mom_p);
   t->SetBranchAddress("mom_pt", &mom_pt, &b_mom_pt);
   t->SetBranchAddress("mom_eta", &mom_eta, &b_mom_eta);
   t->SetBranchAddress("mom_phi", &mom_phi, &b_mom_phi);
   t->SetBranchAddress("mom_vx", &mom_vx, &b_mom_vx);
   t->SetBranchAddress("mom_vy", &mom_vy, &b_mom_vy);
   t->SetBranchAddress("mom_vz", &mom_vz, &b_mom_vz);
   t->SetBranchAddress("mom_beta", &mom_beta, &b_mom_beta);
   t->SetBranchAddress("mom_gama", &mom_gama, &b_mom_gama);
   t->SetBranchAddress("mom_Lxy", &mom_Lxy, &b_mom_Lxy);
   t->SetBranchAddress("mom_Lz", &mom_Lz, &b_mom_Lz);
   t->SetBranchAddress("mom_Lxyz", &mom_Lxyz, &b_mom_Lxyz);
   t->SetBranchAddress("mom_ctau", &mom_ctau, &b_mom_ctau);
   t->SetBranchAddress("mom_dupl", &mom_dupl, &b_mom_dupl);

   // setup histo map
   TH1map h1map;
   TH2map h2map;
   histos(h1map, h2map);
   
   int good_ex = 0;

   // loop over events
   unsigned int nentries = t->GetEntries();
   for (unsigned int entry = 0; entry < nentries; entry++){
     t->GetEntry(entry);

     int test_nconst = 0;
     for (unsigned int gp = 0; gp < ngenpart; gp++){
       if ( (*genpar_stat)[gp] != 1 ) continue;           // keep only final state particles
       if ( std::abs((*genpar_eta)[gp]) > 1.5 ) continue; // keep only particles in barrel 
       if ( (*genpar_pt)[gp] < 1.0 ) continue;            // keep only >1 GeV particles
       if ( (*genpar_match_q1)[gp]!=1) continue;          // keep only q1 particles
       test_nconst++;
     }
     if (test_nconst>20) good_ex += 1;
  
     std::vector<int> same_mom;
     same_mom.resize(4);
     for (unsigned int gm = 0; gm < 4; gm++){
       if ((*mom_id)[gm]==35)
       if ((*mom_dupl)[gm]==1) continue; // only plot mom once
       h1map["LL_beta"]->Fill((*mom_beta)[gm]);
       h1map["LL_cTau"]->Fill((*mom_ctau)[gm]);

     }

     float jet1_mom_beta   = (*mom_beta)[0];
     float jet1_mom_lxyz   = (*mom_Lxyz)[0];
     float jet2_mom_beta   = (*mom_beta)[1];
     float jet2_mom_lxyz   = (*mom_Lxyz)[1];
     float jet3_mom_beta   = (*mom_beta)[2];
     float jet3_mom_lxyz   = (*mom_Lxyz)[2];
     float jet4_mom_beta   = (*mom_beta)[3];
     float jet4_mom_lxyz   = (*mom_Lxyz)[3];

     float jet_const_beta  = 0;
     float jet_const_lxyz  = 0;
     float jet_orig_beta   = 0;
     float jet_orig_lxyz   = 0;

     float jet1_const_dt   = 0;
     float jet1_time_raw   = 0;
     float jet1_nconst     = 0;
     float jet2_const_dt   = 0;
     float jet2_time_raw   = 0;
     float jet2_nconst     = 0;
     float jet3_const_dt   = 0;
     float jet3_time_raw   = 0;
     float jet3_nconst     = 0;
     float jet4_const_dt   = 0;
     float jet4_time_raw   = 0;
     float jet4_nconst     = 0;
     float unmatch_dt      = 0;

     // loop over gen particles

     for (unsigned int gp = 0; gp < ngenpart; gp++){

       if ( (*genpar_stat)[gp] != 1 ) continue;           // keep only final state particles
       if ( std::abs((*genpar_eta)[gp]) > 1.5 ) continue; // keep only particles in barrel 
       if ( (*genpar_pt)[gp] < 1.0 ) continue;            // keep only >1 GeV particles

       // check particle matches jet
       if ( (*genpar_match_q1)[gp]==1 || (*genpar_match_q2)[gp]==1 || (*genpar_match_q3)[gp]==1 || (*genpar_match_q4)[gp]==1 ){ 
         // constituent info
         h1map["const_id"]->Fill((*genpar_id)[gp]);
         h1map["const_pt"]->Fill((*genpar_pt)[gp]);
         h1map["const_eta"]->Fill(std::abs((*genpar_eta)[gp]));
       }
   
       // constituent particle info
       jet_const_beta = 1.0;
       jet_const_lxyz = (*genpar_la)[gp];
       // hypothetical non-LL info
       jet_orig_beta  = 1.0;
       jet_orig_lxyz  = (*genpar_lo)[gp];

       // --- q1 jet
       if ( (*genpar_match_q1)[gp]==1){
         jet1_nconst  += 1.0; // counter for jet constituents
         jet1_const_dt = calcDeltaT(jet1_mom_lxyz,jet1_mom_beta,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta); // constituent deltaT 
         jet1_time_raw += jet1_const_dt;
         h1map["const_t"]->Fill(jet1_const_dt);
         if (good_ex==1) h1map["const_t_ex"]->Fill(jet1_const_dt);
       }// end match to q1 jet

       // --- q2 jet
       if ( (*genpar_match_q2)[gp]==1){
         jet2_nconst  += 1.0; // counter for jet constituents
         jet2_const_dt = calcDeltaT(jet2_mom_lxyz,jet2_mom_beta,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta); // constituent deltaT 
         jet2_time_raw += jet2_const_dt;
         h1map["const_t"]->Fill(jet2_const_dt);
       }// end match to q2 jet

       // --- q3 jet
       if ( (*genpar_match_q3)[gp]==1){
         jet3_nconst  += 1.0; // counter for jet constituents
         jet3_const_dt = calcDeltaT(jet3_mom_lxyz,jet3_mom_beta,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta); // constituent deltaT 
         jet3_time_raw += jet3_const_dt;
         h1map["const_t"]->Fill(jet3_const_dt);
       }// end match to q3 jet

       // --- q4 jet
       if ( (*genpar_match_q4)[gp]==1){
         jet4_nconst  += 1.0; // counter for jet constituents
         jet4_const_dt = calcDeltaT(jet4_mom_lxyz,jet4_mom_beta,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta); // constituent deltaT 
         jet4_time_raw += jet4_const_dt;
         h1map["const_t"]->Fill(jet4_const_dt);
       }// end match to q4 jet

       // --- unmatched particles
       //if ( (*genpar_match_q1)[gp]==0 && (*genpar_match_q2)[gp]==0 && (*genpar_match_q3)[gp]==0 && (*genpar_match_q4)[gp]==0 ){
       if ( (*genpar_match_q0)[gp]==0 ){
         unmatch_dt = calcDeltaT(0.0,1.0,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta);
         h1map["unmatch_t"]->Fill(unmatch_dt);
       }// end unmatched particles

     }// end loop over gen particles

     float jet1_time_avg = calcAvgT(jet1_time_raw,jet1_nconst);
     float jet2_time_avg = calcAvgT(jet2_time_raw,jet2_nconst);
     float jet3_time_avg = calcAvgT(jet3_time_raw,jet3_nconst);
     float jet4_time_avg = calcAvgT(jet4_time_raw,jet4_nconst);

     if (jet1_nconst > 0) h1map["nconst"]->Fill(jet1_nconst);
     if (jet2_nconst > 0) h1map["nconst"]->Fill(jet2_nconst);
     if (jet3_nconst > 0) h1map["nconst"]->Fill(jet3_nconst);
     if (jet4_nconst > 0) h1map["nconst"]->Fill(jet4_nconst);
      
     if (jet1_nconst > 0) h1map["jet_t"]->Fill(jet1_time_avg);
     if (jet2_nconst > 0) h1map["jet_t"]->Fill(jet2_time_avg);
     if (jet3_nconst > 0) h1map["jet_t"]->Fill(jet3_time_avg);
     if (jet4_nconst > 0) h1map["jet_t"]->Fill(jet4_time_avg);

     float max_12 = std::max(jet1_time_avg,jet2_time_avg);
     float max_34 = std::max(jet3_time_avg,jet4_time_avg);
     float max_jet_t = std::max(max_12,max_34);
     h1map["max_jet_t"]->Fill(max_jet_t);

     float dt_jetsX1   = 0;
     float dt_jetsX2   = 0;
     float dt_jet_near = 0;
     float dt_jet_far  = 0;
     h2map["tX1_tX2"]->Fill(dt_jetsX1,dt_jetsX2);
     h2map["tnear_tfar"]->Fill(dt_jet_near,dt_jet_far);

  }// end loop over events
 
  // save plots
  save1Dplots(out,fout,h1map);
  save2Dplots(out,fout,h2map);

}// end run()

float calcDeltaT(const float lx, const float bx, const float la, const float ba, const float lo, const float bo){
  // change in distance
  float dl = (float)lx/(float)bx + (float)la/(float)ba - (float)lo/(float)bo;
  // convert to time by dividing by c (30cm/ns)
  float dt = (float)dl/30.0; 
  return dt;
}// end calcDeltaT

float calcAvgT(const float t, const int n){
  float avgT = -1000;
  if (n!=0) avgT = (float)t/(float)n;
  return avgT;
}// end calcAvgT

