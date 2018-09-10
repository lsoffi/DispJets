#include <iostream>

typedef std::map<TString, TH1F*> TH1map;
typedef std::map<TString, TH2F*> TH2map;

// declare other functions
void run( TString, TString, TFile* );
void histos( TH1map & , TH2map & );
void save1Dplots( TString, TFile* , const TH1map & );
void save2Dplots( TString, TFile* , const TH2map & );
TH1F * MakeTH1FPlot( const TString, const TString, const int, const double, const double, const TString, const TString );
TH2F * MakeTH2FPlot( const TString, const TString, const int, const double, const double, const int, const double, const double, const TString, const TString );

// START
void drawFromTree(){

  TString path = "../samples/";

  std::vector< TString > ctau;
  ctau.push_back("0mm");
  ctau.push_back("1mm");
  ctau.push_back("10mm");
  ctau.push_back("100mm");
  ctau.push_back("1000mm");
  ctau.push_back("10000mm");
  std::vector< TString > file;
  //file.push_back(Form("%sntuple_dispjets.root",path.Data()));
  for (int i = 0; i < ctau.size(); i++){
    file.push_back(Form("%sntuple_xxqqqq_m50_ct%s.root",path.Data(),ctau[i].Data()));
  }
  //file.push_back(Form("%sntuple_xxqqqq_m50_ct0mm.root",path.Data()));
  //file.push_back(Form("%sntuple_xxqqqq_m50_ct1mm.root",path.Data()));
  //file.push_back(Form("%sntuple_xxqqqq_m50_ct10mm.root",path.Data()));
  //file.push_back(Form("%sntuple_xxqqqq_m50_ct100mm.root",path.Data()));
  //file.push_back(Form("%sntuple_xxqqqq_m50_ct1000mm.root",path.Data()));
  //file.push_back(Form("%sntuple_xxqqqq_m50_ct10000mm.root",path.Data()));
  int nsamples = file.size();

  for (int f = 0; f < nsamples; f++){
    TString out  = "~/www/Plots/DispJets/GenLevelPlots/NewSamples_v3/ct"+ctau[f]+"/";
    TFile *fout = TFile::Open(Form("%splots_dispjets.root",out.Data()),"RECREATE");
    run(file[f],out,fout);
  } 

}

void histos( TH1map & map , TH2map & map2){

  // 1D histos
  map["ngj"]		= MakeTH1FPlot("ngj","",20,0,20,"Num. gen jets","");
  map["ngp"]		= MakeTH1FPlot("ngp","",20,0,20,"Num. gen particles","");

  // gen jet properties
  map["gj_pt"]		= MakeTH1FPlot("gj_pt","",20,0,100,"Gen jet p_T","");
  map["gj_e"]		= MakeTH1FPlot("gj_e","",20,0,100,"Gen jet E",""); 
  map["gj_eta"]		= MakeTH1FPlot("gj_eta","",20,-4,4,"Gen jet #eta","");
  map["gj_phi"]		= MakeTH1FPlot("gj_phi","",20,-4,4,"Gen jet #phi","");
  map["gj_nd"]		= MakeTH1FPlot("gj_nd","",10,0,10,"Gen jet num. daughters","");
  map["gj_nc"]		= MakeTH1FPlot("gj_nc","",10,0,10,"Gen jet num. constituents","");
  map["gj_vx"]		= MakeTH1FPlot("gj_vx","",400,-200,200,"Gen jet vertex X","");
  map["gj_vy"]		= MakeTH1FPlot("gj_vy","",400,-200,200,"Gen jet vertex Y","");
  map["gj_vz"]		= MakeTH1FPlot("gj_vz","",400,-200,200,"Gen jet vertex Z","");
  map["gj0_c_id"]	= MakeTH1FPlot("gj0_c_id","",400,-200,200,"Gen jet 0 const. ID","");
  map["gj0_c_pt"]	= MakeTH1FPlot("gj0_c_pt","",300,0,300,"Gen jet 0 const. p_T","");
  map["gj0_c_vx"]	= MakeTH1FPlot("gj0_c_vx","",400,-200,200,"Gen jet 0 const. vertex X","");
  map["gj0_c_vy"]	= MakeTH1FPlot("gj0_c_vy","",400,-200,200,"Gen jet 0 const. vertex Y","");
  map["gj0_c_vz"]	= MakeTH1FPlot("gj0_c_vz","",400,-200,200,"Gen jet 0 const. vertex Z","");
  map["gj1_c_id"]	= MakeTH1FPlot("gj1_c_id","",400,-200,200,"Gen jet 1 const. ID","");
  map["gj1_c_pt"]	= MakeTH1FPlot("gj1_c_pt","",300,0,300,"Gen jet 1 const. p_T","");
  map["gj1_c_vx"]	= MakeTH1FPlot("gj1_c_vx","",400,-200,200,"Gen jet 1 const. vertex X","");
  map["gj1_c_vy"]	= MakeTH1FPlot("gj1_c_vy","",400,-200,200,"Gen jet 1 const. vertex Y","");
  map["gj1_c_vz"]	= MakeTH1FPlot("gj1_c_vz","",400,-200,200,"Gen jet 1 const. vertex Z","");
  map["gj2_c_id"]	= MakeTH1FPlot("gj2_c_id","",400,-200,200,"Gen jet 2 const. ID","");
  map["gj2_c_pt"]	= MakeTH1FPlot("gj2_c_pt","",300,0,300,"Gen jet 2 const. p_T","");
  map["gj2_c_vx"]	= MakeTH1FPlot("gj2_c_vx","",400,-200,200,"Gen jet 2 const. vertex X","");
  map["gj2_c_vy"]	= MakeTH1FPlot("gj2_c_vy","",400,-200,200,"Gen jet 2 const. vertex Y","");
  map["gj2_c_vz"]	= MakeTH1FPlot("gj2_c_vz","",400,-200,200,"Gen jet 2 const. vertex Z","");
  map["gj3_c_id"]	= MakeTH1FPlot("gj3_c_id","",400,-200,200,"Gen jet 3 const. ID","");
  map["gj3_c_pt"]	= MakeTH1FPlot("gj3_c_pt","",300,0,300,"Gen jet 3 const. p_T","");
  map["gj3_c_vx"]	= MakeTH1FPlot("gj3_c_vx","",400,-200,200,"Gen jet 3 const. vertex X","");
  map["gj3_c_vy"]	= MakeTH1FPlot("gj3_c_vy","",400,-200,200,"Gen jet 3 const. vertex Y","");
  map["gj3_c_vz"]	= MakeTH1FPlot("gj3_c_vz","",400,-200,200,"Gen jet 3 const. vertex Z","");

  // gen particle properties
  map["gp_id"]		= MakeTH1FPlot("gp_id","",200,-100,100,"Gen particle ID","");
  map["gp_pt"]		= MakeTH1FPlot("gp_pt","",20,0,100,"Gen particle p_T","");
  map["gp_e"]		= MakeTH1FPlot("gp_e","",20,0,100,"Gen particle E",""); 
  map["gp_eta"]		= MakeTH1FPlot("gp_eta","",20,-4,4,"Gen particle #eta","");
  map["gp_phi"]		= MakeTH1FPlot("gp_phi","",20,-4,4,"Gen particle #phi","");
  map["gp_vx"]		= MakeTH1FPlot("gp_vx","",400,-200,200,"Gen particle vertex X","");
  map["gp_vy"]		= MakeTH1FPlot("gp_vy","",400,-200,200,"Gen particle vertex Y","");
  map["gp_vz"]		= MakeTH1FPlot("gp_vz","",400,-200,200,"Gen particle vertex Z","");
  map["gp_Lxy"]		= MakeTH1FPlot("gp_Lxy","",100,0,100,"Gen particle vertex XY dist","");
  map["gp_Lxyz"]	= MakeTH1FPlot("gp_Lxyz","",100,0,100,"Gen particle vertex XYZ dist","");

  // gen particle mother properties
  map["gm_id"]		= MakeTH1FPlot("gm_id","",200,-100,100,"Gen mom ID","");
  map["gm_stat"]	= MakeTH1FPlot("gm_stat","",200,-100,100,"Gen mom status","");
  map["gm_pt"]		= MakeTH1FPlot("gm_pt","",20,0,100,"Gen mom p_T","");
  map["gm_e"]		= MakeTH1FPlot("gm_e","",20,0,100,"Gen mom E","");
  map["gm_m"]		= MakeTH1FPlot("gm_m","",100,0,500,"Gen mom mass",""); 
  map["gm_eta"]		= MakeTH1FPlot("gm_eta","",20,-4,4,"Gen mom #eta","");
  map["gm_phi"]		= MakeTH1FPlot("gm_phi","",20,-4,4,"Gen mom #phi","");
  map["gm_vx"]		= MakeTH1FPlot("gm_vx","",400,-200,200,"Gen mom vertex X","");
  map["gm_vy"]		= MakeTH1FPlot("gm_vy","",400,-200,200,"Gen mom vertex Y","");
  map["gm_vz"]		= MakeTH1FPlot("gm_vz","",400,-200,200,"Gen mom vertex Z","");
  map["gm_Lxy"]		= MakeTH1FPlot("gm_Lxy","",100,0,100,"Gen mom vertex XY dist","");
  map["gm_Lz"]		= MakeTH1FPlot("gm_Lz","",100,0,100,"Gen mom vertex Z dist","");
  map["gm_Lxyz"]	= MakeTH1FPlot("gm_Lxyz","",100,0,100,"Gen mom vertex XYZ dist","");
  map["gm_ctau"]	= MakeTH1FPlot("gm_ctau","",100,0,500,"Gen mom c#tau","");

  // 2D histos
  map2["gp_vx_id"]	= MakeTH2FPlot("gp_vx_id","",400,-200,200,400,-200,200,"Gen particle ID","Gen particle vertex X");
  map2["gp_vy_id"]	= MakeTH2FPlot("gp_vy_id","",400,-200,200,400,-200,200,"Gen particle ID","Gen particle vertex Y");
  map2["gp_vz_id"]	= MakeTH2FPlot("gp_vz_id","",400,-200,200,400,-200,200,"Gen particle ID","Gen particle vertex Z");

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
  vector<int>     *genjet_i;

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
  vector<float>   *genpar_pt;
  vector<float>   *genpar_eta;
  vector<float>   *genpar_phi;
  vector<float>   *genpar_vx;
  vector<float>   *genpar_vy;
  vector<float>   *genpar_vz;
  vector<float>   *genpar_Lxy;
  vector<float>   *genpar_Lxyz;
  vector<int>     *mom_id;
  vector<int>     *mom_stat;
  vector<float>   *mom_e;
  vector<float>   *mom_m;
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

   // Set object pointer
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
   genpar_pt = 0;
   genpar_eta = 0;
   genpar_phi = 0;
   genpar_vx = 0;
   genpar_vy = 0;
   genpar_vz = 0;
   genpar_Lxy = 0;
   genpar_Lxyz = 0;
   mom_id = 0;
   mom_stat = 0;
   mom_e = 0;
   mom_m = 0;
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

  // list of branches
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
  TBranch        *b_genpar_pt;   //!
  TBranch        *b_genpar_eta;   //!
  TBranch        *b_genpar_phi;   //!
  TBranch        *b_genpar_vx;   //!
  TBranch        *b_genpar_vy;   //!
  TBranch        *b_genpar_vz;   //!
  TBranch        *b_genpar_Lxy;   //!
  TBranch        *b_genpar_Lxyz;   //!
  TBranch        *b_mom_id;   //!
  TBranch        *b_mom_stat;   //!
  TBranch        *b_mom_e;   //!
  TBranch        *b_mom_m;   //!
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

  // setup branches
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
  t->SetBranchAddress("genpar_pt", &genpar_pt, &b_genpar_pt);
  t->SetBranchAddress("genpar_eta", &genpar_eta, &b_genpar_eta);
  t->SetBranchAddress("genpar_phi", &genpar_phi, &b_genpar_phi);
  t->SetBranchAddress("genpar_vx", &genpar_vx, &b_genpar_vx);
  t->SetBranchAddress("genpar_vy", &genpar_vy, &b_genpar_vy);
  t->SetBranchAddress("genpar_vz", &genpar_vz, &b_genpar_vz);
  t->SetBranchAddress("genpar_Lxy", &genpar_Lxy, &b_genpar_Lxy);
  t->SetBranchAddress("genpar_Lxyz", &genpar_Lxyz, &b_genpar_Lxyz);
  t->SetBranchAddress("mom_id", &mom_id, &b_mom_id);
  t->SetBranchAddress("mom_stat", &mom_stat, &b_mom_stat);
  t->SetBranchAddress("mom_e", &mom_e, &b_mom_e);
  t->SetBranchAddress("mom_m", &mom_m, &b_mom_m);
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

    h1map["ngj"]->Fill(ngenjets);
    h1map["ngp"]->Fill(ngenpart); 

    // loop over gen jets
    for (unsigned int gj = 0; gj < ngenjets; gj++){
      h1map["gj_pt"]->Fill((*genjet_pt)[gj]);
      h1map["gj_e"]->Fill((*genjet_e)[gj]);
      h1map["gj_eta"]->Fill((*genjet_eta)[gj]);
      h1map["gj_phi"]->Fill((*genjet_phi)[gj]);
      h1map["gj_nd"]->Fill((*genjet_ndaug)[gj]);
      h1map["gj_nc"]->Fill((*genjet_nconst)[gj]);
      h1map["gj_vx"]->Fill((*genjet_vx)[gj]); 
      h1map["gj_vy"]->Fill((*genjet_vy)[gj]); 
      h1map["gj_vz"]->Fill((*genjet_vz)[gj]); 
    }// end loop over genjets

    // loop over particles in gen jet0
    for (unsigned int i = 0; i < (*genjet0_const_id).size(); i++){
      if ((*genjet0_const_st)[i] != 1) continue;
      h1map["gj0_c_id"]->Fill((*genjet0_const_id)[i]);
      h1map["gj0_c_pt"]->Fill((*genjet0_const_pt)[i]);
      h1map["gj0_c_vx"]->Fill((*genjet0_const_vx)[i]);
      h1map["gj0_c_vy"]->Fill((*genjet0_const_vy)[i]);
      h1map["gj0_c_vz"]->Fill((*genjet0_const_vz)[i]);
    }
    // loop over particles in gen jet1
    for (unsigned int i = 0; i < (*genjet1_const_id).size(); i++){
      if ((*genjet1_const_st)[i] != 1) continue;
      h1map["gj1_c_id"]->Fill((*genjet1_const_id)[i]);
      h1map["gj1_c_pt"]->Fill((*genjet1_const_pt)[i]);
      h1map["gj1_c_vx"]->Fill((*genjet1_const_vx)[i]);
      h1map["gj1_c_vy"]->Fill((*genjet1_const_vy)[i]);
      h1map["gj1_c_vz"]->Fill((*genjet1_const_vz)[i]);
    }
    // loop over particles in gen jet2
    for (unsigned int i = 0; i < (*genjet2_const_id).size(); i++){
      if ((*genjet2_const_st)[i] != 1) continue;
      h1map["gj2_c_id"]->Fill((*genjet2_const_id)[i]);
      h1map["gj2_c_pt"]->Fill((*genjet2_const_pt)[i]);
      h1map["gj2_c_vx"]->Fill((*genjet2_const_vx)[i]);
      h1map["gj2_c_vy"]->Fill((*genjet2_const_vy)[i]);
      h1map["gj2_c_vz"]->Fill((*genjet2_const_vz)[i]);
    }
    // loop over particles in gen jet3
    for (unsigned int i = 0; i < (*genjet3_const_id).size(); i++){
      if ((*genjet3_const_st)[i] != 1) continue;
      h1map["gj3_c_id"]->Fill((*genjet3_const_id)[i]);
      h1map["gj3_c_pt"]->Fill((*genjet3_const_pt)[i]);
      h1map["gj3_c_vx"]->Fill((*genjet3_const_vx)[i]);
      h1map["gj3_c_vy"]->Fill((*genjet3_const_vy)[i]);
      h1map["gj3_c_vz"]->Fill((*genjet3_const_vz)[i]);
    }

    // loop over gen particles
    for (unsigned int gp = 0; gp < ngenpart; gp++){
      h2map["gp_vx_id"]->Fill((*genpar_id)[gp],(*genpar_vx)[gp]);
      h2map["gp_vy_id"]->Fill((*genpar_id)[gp],(*genpar_vy)[gp]);
      h2map["gp_vz_id"]->Fill((*genpar_id)[gp],(*genpar_vz)[gp]);
      // gen particles info
      h1map["gp_id"]->Fill((*genpar_id)[gp]);
      h1map["gp_pt"]->Fill((*genpar_pt)[gp]);
      h1map["gp_eta"]->Fill((*genpar_eta)[gp]);
      h1map["gp_phi"]->Fill((*genpar_phi)[gp]);
      h1map["gp_vx"]->Fill((*genpar_vx)[gp]);
      h1map["gp_vy"]->Fill((*genpar_vy)[gp]);
      h1map["gp_vz"]->Fill((*genpar_vz)[gp]);
      h1map["gp_Lxy"]->Fill((*genpar_Lxy)[gp]);
      h1map["gp_Lxyz"]->Fill((*genpar_Lxyz)[gp]);
      if ((*mom_id)[gp] == 0 || (*mom_stat)[gp] == 0) continue;
      //if (entry==0) std::cout << "ID: " << (*genpar_id)[gp] << " mom ID: " << (*mom_id)[gp] << std::endl;
      // mother particle info
      h1map["gm_id"]->Fill((*mom_id)[gp]);
      h1map["gm_stat"]->Fill((*mom_stat)[gp]);
      h1map["gm_e"]->Fill((*mom_e)[gp]);
      h1map["gm_m"]->Fill((*mom_m)[gp]);
      h1map["gm_pt"]->Fill((*mom_pt)[gp]);
      h1map["gm_eta"]->Fill((*mom_eta)[gp]);
      h1map["gm_phi"]->Fill((*mom_phi)[gp]);
      h1map["gm_vx"]->Fill((*mom_vx)[gp]);
      h1map["gm_vy"]->Fill((*mom_vy)[gp]);
      h1map["gm_vz"]->Fill((*mom_vz)[gp]);
      h1map["gm_Lxy"]->Fill((*mom_Lxy)[gp]);
      h1map["gm_Lz"]->Fill((*mom_Lz)[gp]);
      h1map["gm_Lxyz"]->Fill((*mom_Lxyz)[gp]);
      h1map["gm_ctau"]->Fill((*mom_ctau)[gp]);
    }// end loop over gen particles

  }// end loop over events
 
  // save plots
  save1Dplots(out,fout,h1map);
  save2Dplots(out,fout,h2map);
 
} 

