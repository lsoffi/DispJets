#include <iostream>

typedef std::map<TString, TH1F*>	TH1map;
typedef TH1map::iterator		TH1mapIter;


// declare other functions
void run( TString, TString, TString, TFile* );
void histos( TH1map & );
void saveplots( TString, TFile* , TH1map & );
TH1F * MakeTH1FPlot(const TString, const TString, const int, const double, const double, const TString, const TString);


// START
void drawFromTree(){

  TString path = "../";
  TString out  = "~/www/Plots/DispJets/GenLevelPlots/";
  std::vector< TString > file;
  file.push_back("ntuple_dispjets.root");
  int nsamples = file.size();

  TFile *fout = TFile::Open(Form("%splots_dispjets.root",out.Data()),"RECREATE");
 
  for (int f = 0; f < nsamples; f++){
    run(path,file[f],out,fout);
  } 

}

void histos( TH1map & map ){

  map["ngj"]			= MakeTH1FPlot("ngj","",20,0,20,"Num. gen jets","");
  map["ngp"]			= MakeTH1FPlot("ngp","",20,0,20,"Num. gen particles","");
  map["gj_vx"]			= MakeTH1FPlot("gj_vx","",100,-50,50,"Gen jet vertex X","");
  map["gj_vy"]			= MakeTH1FPlot("gj_vy","",100,-50,50,"Gen jet vertex Y","");
  map["gj_vz"]			= MakeTH1FPlot("gj_vz","",100,-50,50,"Gen jet vertex Z","");
  map["gp_vx"]			= MakeTH1FPlot("gp_vx","",100,-50,50,"Gen particle vertex X","");
  map["gp_vy"]			= MakeTH1FPlot("gp_vy","",100,-50,50,"Gen particle vertex Y","");
  map["gp_vz"]			= MakeTH1FPlot("gp_vz","",100,-50,50,"Gen particle vertex Z","");
  map["gm_vx"]			= MakeTH1FPlot("gm_vx","",100,-50,50,"Gen mom vertex X","");
  map["gm_vy"]			= MakeTH1FPlot("gm_vy","",100,-50,50,"Gen mom vertex Y","");
  map["gm_vz"]			= MakeTH1FPlot("gm_vz","",100,-50,50,"Gen mom vertex Z","");

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

void saveplots(TString odir, TFile* fout, TH1map & map){
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

}// end saveplots

void run(TString path, TString file, TString out, TFile* fout){

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
  TH1map hmap;
  histos(hmap);

  // loop over events
  unsigned int nentries = t->GetEntries();
  for (unsigned int entry = 0; entry < nentries; entry++){
    t->GetEntry(entry);

    hmap["ngj"]->Fill(ngenjets);
    hmap["ngp"]->Fill(ngenpart); 

    // loop over gen jets
    for (unsigned int gj = 0; gj < ngenjets; gj++){
      hmap["gj_vx"]->Fill((*genjet_vx)[gj]); 
      hmap["gj_vy"]->Fill((*genjet_vy)[gj]); 
      hmap["gj_vz"]->Fill((*genjet_vz)[gj]); 
    }// end loop over genjets

    // loop over gen particles
    for (unsigned int gp = 0; gp < ngenpart; gp++){
      // gen particles info
      hmap["gp_vx"]->Fill((*genpar_vx)[gp]);
      hmap["gp_vy"]->Fill((*genpar_vy)[gp]);
      hmap["gp_vz"]->Fill((*genpar_vz)[gp]);
      // mother particle info
      hmap["gm_vx"]->Fill((*mom_vx)[gp]);
      hmap["gm_vy"]->Fill((*mom_vy)[gp]);
      hmap["gm_vz"]->Fill((*mom_vz)[gp]);
    }// end loop over gen particles

  }// end loop over events
 
  // save plots
  saveplots(out,fout,hmap);
 
} 

