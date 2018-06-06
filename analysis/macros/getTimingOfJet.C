#include <iostream>

// declare types
typedef std::map<TString, TH1F*> TH1map;
typedef std::map<TString, TH2F*> TH2map;

// declare other functions
void getJetTime( );
void run( TString, TString, TFile* );
void histos( TH1map & , TH2map & );
void save1Dplots( TString, TFile* , const TH1map & );
void save2Dplots( TString, TFile* , const TH2map & );
TH1F * MakeTH1FPlot( const TString, const TString, const int, const double, const double, const TString, const TString );
TH2F * MakeTH2FPlot( const TString, const TString, const int, const double, const double, const int, const double, const double, const TString, const TString );


// ---------------------------------- 
// ----------- START ----------------
// ---------------------------------- 
void getTimingOfJet(){

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
  map["LL_cTau"]	= MakeTH1FPlot("LL_cTau","",50,0,50,"LL particle c#tau","");
  map["const_t_ev1"]	= MakeTH1FPlot("const_t_ev1","",200,0,20,"Jet constituent time [ns] Ev.1","");
  map["nconst"]		= MakeTH1FPlot("nconst","",50,0,50,"Num. jet constituents","");
  map["const_t"]	= MakeTH1FPlot("const_t","",200,0,20,"Jet constituent time [ns]",""); 
  map["const_id"]	= MakeTH1FPlot("const_id","",400,-200,200,"Jet constituent pdgId","");
  map["const_pt"]	= MakeTH1FPlot("const_pt","",300,0,300,"Jet constituent p_T","");
  map["const_eta"]	= MakeTH1FPlot("const_eta","",30,0,3.0,"Jet constituent |#eta|","");
  map["jet_t"]		= MakeTH1FPlot("jet_t","",200,0,20,"Jet time [ns]",""); 

}



TH1F * MakeTH1FPlot(const TString hname, const TString htitle, const int nbins, const double xlow, const double xhigh, const TString xtitle, const TString ytitle){
  TString ytitleNew;
  Float_t binwidth = (xhigh-xlow)/nbins;
  if (ytitle=="") ytitleNew = Form("Events / %2.1f GeV",binwidth);
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
    c->SaveAs(Form("%s%s.png",odir.Data(),name.Data()));
    c->SaveAs(Form("%s%s.pdf",odir.Data(),name.Data()));

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
   vector<int>     *genpar_match_q1;
   vector<int>     *genpar_match_q2;
   vector<int>     *genpar_match_q3;
   vector<int>     *genpar_match_q4;
   vector<float>   *genpar_lo;
   vector<float>   *genpar_la;
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
   TBranch        *b_genpar_match_q1;   //!
   TBranch        *b_genpar_match_q2;   //!
   TBranch        *b_genpar_match_q3;   //!
   TBranch        *b_genpar_match_q4;   //!
   TBranch        *b_genpar_lo;   //!
   TBranch        *b_genpar_la;   //!
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
   t->SetBranchAddress("genpar_match_q1", &genpar_match_q1, &b_genpar_match_q1);
   t->SetBranchAddress("genpar_match_q2", &genpar_match_q2, &b_genpar_match_q2);
   t->SetBranchAddress("genpar_match_q3", &genpar_match_q3, &b_genpar_match_q3);
   t->SetBranchAddress("genpar_match_q4", &genpar_match_q4, &b_genpar_match_q4);
   t->SetBranchAddress("genpar_lo", &genpar_lo, &b_genpar_lo);
   t->SetBranchAddress("genpar_la", &genpar_la, &b_genpar_la);
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

   // setup histo map
   TH1map h1map;
   TH2map h2map;
   histos(h1map, h2map);
   
   // loop over events
   unsigned int nentries = t->GetEntries();
   for (unsigned int entry = 0; entry < nentries; entry++){
     t->GetEntry(entry);

     int nmothers = (*mom_beta).size();

     std::cout << "nmom " << nmothers << std::endl;

     for (unsigned int gm = 0; gm < nmothers; gm++){
       h1map["LL_beta"]->Fill((*mom_beta)[gm]);
       h1map["LL_cTau"]->Fill((*mom_ctau)[gm]);
     }

     float jet1_mom_beta   = (*mom_beta)[1];
     float jet1_mom_lxyz   = (*mom_Lxyz)[1];

     float jet1_const_beta = 0;
     float jet1_const_lxyz = 0;
     float jet1_orig_beta  = 0;
     float jet1_orig_lxyz  = 0;

     float jet1_const_dl   = 0;
     float jet1_const_dt   = 0;
     float jet1_time_raw   = 0;
     float jet1_nconst     = 0;
     float jet1_time_avg   = 0;

     // loop over gen particles
     for (unsigned int gp = 0; gp < ngenpart; gp++){

       if ( (*genpar_stat)[gp] != 1 ) continue; // keep only final state particles
       
       if ( (*genpar_match_q1)[gp] == 1){ // check particle matches q1 jet

         // constituent info
         h1map["const_id"]->Fill((*genpar_id)[gp]);
         h1map["const_pt"]->Fill((*genpar_pt)[gp]);
         h1map["const_eta"]->Fill(std::abs((*genpar_eta)[gp]));

         jet1_nconst+=1.0; // counter for jet constituents
         
         // const particle
         jet1_const_beta = 1.0;
         jet1_const_lxyz = (*genpar_la)[gp];

         // not-LL 
         jet1_orig_beta  = 1.0;
         jet1_orig_lxyz  = (*genpar_lo)[gp];

         // delta T for const
         jet1_const_dl  = (float)jet1_mom_lxyz/(float)jet1_mom_beta + (float)jet1_const_lxyz/(float)jet1_const_beta - (float)jet1_orig_lxyz/(float)jet1_orig_beta;
         jet1_const_dt  = (float)jet1_const_dl/30.0; // divide by c (30cm/ns)
         h1map["const_t"]->Fill(jet1_const_dt);
         if (entry==0) h1map["const_t_ev1"]->Fill(jet1_const_dt);
         // total delta T for jet
         jet1_time_raw   += jet1_const_dt;


       }// end match to q1 jet

     }// end loop over gen particles

     jet1_time_avg = (float)jet1_time_raw/(float)jet1_nconst; 

     h1map["nconst"]->Fill(jet1_nconst); 
     h1map["jet_t"]->Fill(jet1_time_avg);

  }// end loop over events
 
  // save plots
  save1Dplots(out,fout,h1map);
  save2Dplots(out,fout,h2map);

}// end run()
