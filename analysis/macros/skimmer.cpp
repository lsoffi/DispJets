#include "skimmer.hh"

skimmer::skimmer(TString path, TString sample)
{

  // input file
  std::cout << "Running sample: " << sample << std::endl;
  infile = TFile::Open(Form("%s%s.root",path.Data(),sample.Data()));
  if (!infile){ std::cout << "Sample " << sample << " does not exist! Exiting... " << std::endl; return; }

  // input tree
  intree = (TTree*)infile->Get("dispjets/tree");
  if (!intree){ std::cout << "Tree does not exist! Exiting... " << std::endl; return; }

  // output file
  outfile = TFile::Open(Form("%s%s_skim.root",path.Data(),sample.Data()),"RECREATE");

  // output tree
  outdir  = outfile->mkdir("dispjets"); outdir->cd(); 
  outtree = (TTree*)infile->Get("dispjets/tree");
  outtree = intree->CloneTree(0); // clone original tree

}// end skimmer

skimmer::~skimmer()
{

  // save and delete
  std::cout << "Finishing up" << std::endl;
  outfile->cd();
  outdir->cd();
  outtree->Write();
  delete intree;
  delete infile;
  delete outtree;
  delete outfile; 

}// end ~skimmer


void skimmer::run()
{

  // variables needed to do cuts
  int		ngenpart;					TBranch *b_ngenpart;
  vector<int>	*genpar_stat;		genpar_stat = 0;	TBranch *b_genpar_stat;
  vector<float>	*genpar_eta;		genpar_eta = 0;		TBranch *b_genpar_eta;
  vector<float>	*genpar_pt;		genpar_pt = 0;		TBranch *b_genpar_pt;
  vector<float>	*genpar_la;		genpar_la = 0;		TBranch *b_genpar_la;
  vector<float>	*genpar_lo;		genpar_lo = 0;		TBranch *b_genpar_lo;
  vector<float>	*genpar_beta;		genpar_beta = 0;	TBranch	*b_genpar_beta; 
  vector<float>	*mom_beta;		mom_beta = 0;		TBranch *b_mom_beta;
  vector<float> *mom_Lxyz;		mom_Lxyz = 0;		TBranch *b_mom_Lxyz;
  vector<int>	*genpar_match_q0;	genpar_match_q0 = 0;	TBranch *b_genpar_match_q0;
  vector<int>	*genpar_match_q1;	genpar_match_q1 = 0;	TBranch *b_genpar_match_q1;
  vector<int>	*genpar_match_q2;	genpar_match_q2 = 0;	TBranch *b_genpar_match_q2;
  vector<int>	*genpar_match_q3;	genpar_match_q3 = 0;	TBranch *b_genpar_match_q3;
  vector<int>	*genpar_match_q4;	genpar_match_q4 = 0;	TBranch *b_genpar_match_q4;
  
  intree->SetBranchAddress("ngenpart", 		&ngenpart, 		&b_ngenpart);
  intree->SetBranchAddress("genpar_stat", 	&genpar_stat, 		&b_genpar_stat);
  intree->SetBranchAddress("genpar_eta", 	&genpar_eta, 		&b_genpar_eta);
  intree->SetBranchAddress("genpar_pt", 	&genpar_pt,		&b_genpar_pt);
  intree->SetBranchAddress("genpar_la", 	&genpar_la,		&b_genpar_la);
  intree->SetBranchAddress("genpar_lo", 	&genpar_lo, 		&b_genpar_lo);
  intree->SetBranchAddress("genpar_beta", 	&genpar_beta, 		&b_genpar_beta);
  intree->SetBranchAddress("mom_beta", 		&mom_beta, 		&b_mom_beta);
  intree->SetBranchAddress("mom_Lxyz", 		&mom_Lxyz, 		&b_mom_Lxyz);
  intree->SetBranchAddress("genpar_match_q0",	&genpar_match_q0,	&b_genpar_match_q0); 
  intree->SetBranchAddress("genpar_match_q1",	&genpar_match_q1,	&b_genpar_match_q1); 
  intree->SetBranchAddress("genpar_match_q2",	&genpar_match_q2,	&b_genpar_match_q2); 
  intree->SetBranchAddress("genpar_match_q3",	&genpar_match_q3,	&b_genpar_match_q3); 
  intree->SetBranchAddress("genpar_match_q4",	&genpar_match_q4,	&b_genpar_match_q4); 

  // variables to add to the tree
  vector<int>   jet_nconst;		outtree->Branch("jet_nconst",		&jet_nconst);
  vector<float> jet_avg_t;		outtree->Branch("jet_avg_t",		&jet_avg_t);
  vector<float> jet_smear_30_t;		outtree->Branch("jet_smear_30_t",	&jet_smear_30_t);
  vector<float> jet_smear_180_t;	outtree->Branch("jet_smear_180_t",	&jet_smear_180_t);
  vector<float> jet_pt;			outtree->Branch("jet_pt",		&jet_pt);

  // loop over events
  for (unsigned int entry = 0; entry < intree->GetEntries(); entry++){
    intree->GetEntry(entry);

    // mother info
    float jet1_mom_beta = (*mom_beta)[0];
    float jet1_mom_lxyz = (*mom_Lxyz)[0];
    float jet2_mom_beta = (*mom_beta)[1];
    float jet2_mom_lxyz = (*mom_Lxyz)[1];
    float jet3_mom_beta = (*mom_beta)[2];
    float jet3_mom_lxyz = (*mom_Lxyz)[2];
    float jet4_mom_beta = (*mom_beta)[3];
    float jet4_mom_lxyz = (*mom_Lxyz)[3];

    vector<int>   tmp_jet_nconst;
    vector<float> tmp_jet_pt;	
    vector<float> tmp_jet_avg_t;	
    vector<float> tmp_jet_smear_30_t;	
    vector<float> tmp_jet_smear_180_t;

    int   jet1_nconst         = 0;
    float jet1_pt             = 0;
    float jet1_const_dt       = 0; 
    float jet1_time_raw       = 0;
    float jet1_time_smear_30  = 0; 
    float jet1_time_smear_180 = 0;
    int   jet2_nconst         = 0; 
    float jet2_pt             = 0;
    float jet2_const_dt       = 0; 
    float jet2_time_raw       = 0;
    float jet2_time_smear_30  = 0; 
    float jet2_time_smear_180 = 0;
    int   jet3_nconst         = 0; 
    float jet3_pt             = 0;
    float jet3_const_dt       = 0; 
    float jet3_time_raw       = 0;
    float jet3_time_smear_30  = 0; 
    float jet3_time_smear_180 = 0;
    int   jet4_nconst         = 0; 
    float jet4_pt             = 0;
    float jet4_const_dt       = 0; 
    float jet4_time_raw       = 0;
    float jet4_time_smear_30  = 0; 
    float jet4_time_smear_180 = 0;

    float smear_value_30      = 0; 
    float smear_value_180     = 0;
    float jet_const_beta      = 0; 
    float jet_const_lxyz      = 0;
    float jet_orig_beta       = 0; 
    float jet_orig_lxyz       = 0;  
 

    // loop over gen particles
    for (int gp = 0; gp < ngenpart; gp++){

      if ((*genpar_stat)[gp] != 1)            continue;
      if (std::fabs((*genpar_eta)[gp]) > 1.5) continue;
      if ((*genpar_pt)[gp] < 1.0)             continue;
  
      smear_value_30  = smearVal(0.03); // res: 30 ps -> 0.03 ns
      smear_value_180 = smearVal(0.18); // res: 180 ps -> 0.18 ns
 
      // constituent particle info 
      jet_const_beta = 1.0;
      jet_const_lxyz = (*genpar_la)[gp];
      // hypothetical non-LL info
      jet_orig_beta  = 1.0; // (*genpar_beta)[gp]
      jet_orig_lxyz  = (*genpar_lo)[gp];    
      
      //--- q1 jet
      if ( (*genpar_match_q1)[gp]==1 ){
        jet1_const_dt        = calcDeltaT(jet1_mom_lxyz,jet1_mom_beta,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta);
        jet1_nconst         += 1.0; 
        jet1_time_raw       += jet1_const_dt;
        jet1_time_smear_30  += jet1_const_dt*smear_value_30;
        jet1_time_smear_180 += jet1_const_dt*smear_value_180;
        jet1_pt             += (*genpar_pt)[gp]; 
      }// end match to q1 jet

      // --- q2 jet
      if ( (*genpar_match_q2)[gp]==1){
        jet2_const_dt        = calcDeltaT(jet2_mom_lxyz,jet2_mom_beta,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta);
        jet2_nconst         += 1.0;
        jet2_time_raw       += jet2_const_dt;
        jet2_time_smear_30  += jet2_const_dt*smear_value_30;
        jet2_time_smear_180 += jet2_const_dt*smear_value_180;
        jet2_pt             += (*genpar_pt)[gp]; 
      }// end match to q2 jet

      // --- q3 jet
      if ( (*genpar_match_q3)[gp]==1){
        jet3_const_dt        = calcDeltaT(jet3_mom_lxyz,jet3_mom_beta,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta);
        jet3_nconst         += 1.0;
        jet3_time_raw       += jet3_const_dt;
        jet3_time_smear_30  += jet3_const_dt*smear_value_30;
        jet3_time_smear_180 += jet3_const_dt*smear_value_180;
        jet3_pt             += (*genpar_pt)[gp]; 
      }// end match to q3 jet

      // --- q4 jet
      if ( (*genpar_match_q4)[gp]==1){
        jet4_const_dt        = calcDeltaT(jet4_mom_lxyz,jet4_mom_beta,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta);
        jet4_nconst         += 1.0;
        jet4_time_raw       += jet4_const_dt;
        jet4_time_smear_30  += jet4_const_dt*smear_value_30;
        jet4_time_smear_180 += jet4_const_dt*smear_value_180;
        jet4_pt             += (*genpar_pt)[gp]; 
      }// end match to q4 jet

    }// end loop over gen particles

    // number of constituents
    tmp_jet_nconst.push_back(jet1_nconst);
    tmp_jet_nconst.push_back(jet2_nconst);
    tmp_jet_nconst.push_back(jet3_nconst);
    tmp_jet_nconst.push_back(jet4_nconst);
    // jet pt
    tmp_jet_pt.push_back(jet1_pt);
    tmp_jet_pt.push_back(jet2_pt);
    tmp_jet_pt.push_back(jet3_pt);
    tmp_jet_pt.push_back(jet4_pt);
    // calculate jet time
    tmp_jet_avg_t.push_back(calcAvgT(jet1_time_raw,jet1_nconst));
    tmp_jet_avg_t.push_back(calcAvgT(jet2_time_raw,jet2_nconst));
    tmp_jet_avg_t.push_back(calcAvgT(jet3_time_raw,jet3_nconst));
    tmp_jet_avg_t.push_back(calcAvgT(jet4_time_raw,jet4_nconst));
    // smeared jet time (res 30 ps)
    tmp_jet_smear_30_t.push_back(calcAvgT(jet1_time_smear_30,jet1_nconst));
    tmp_jet_smear_30_t.push_back(calcAvgT(jet2_time_smear_30,jet2_nconst));
    tmp_jet_smear_30_t.push_back(calcAvgT(jet3_time_smear_30,jet3_nconst));
    tmp_jet_smear_30_t.push_back(calcAvgT(jet4_time_smear_30,jet4_nconst));
    // smeared jet time (res 180 ps)
    tmp_jet_smear_180_t.push_back(calcAvgT(jet1_time_smear_180,jet1_nconst));
    tmp_jet_smear_180_t.push_back(calcAvgT(jet2_time_smear_180,jet2_nconst));
    tmp_jet_smear_180_t.push_back(calcAvgT(jet3_time_smear_180,jet3_nconst));
    tmp_jet_smear_180_t.push_back(calcAvgT(jet4_time_smear_180,jet4_nconst));

    // clear new tree vectors 
    jet_nconst.clear();
    jet_pt.clear();
    jet_avg_t.clear();
    jet_smear_30_t.clear();
    jet_smear_180_t.clear();  

    for (int i = 0; i < 4; i++){
      jet_nconst.push_back(tmp_jet_nconst[i]);
      jet_pt.push_back(tmp_jet_pt[i]);
      jet_avg_t.push_back(tmp_jet_avg_t[i]);
      jet_smear_30_t.push_back(tmp_jet_smear_30_t[i]);
      jet_smear_180_t.push_back(tmp_jet_smear_180_t[i]);
    }

    // fill output tree
    outtree->Fill();

  }// end loop over events

}// end run


float skimmer::calcDeltaT(const float lx, const float bx, const float la, const float ba, const float lo, const float bo){
  float dl = (float)lx/(float)bx + (float)la/(float)ba - (float)lo/(float)bo; // change in distance
  float dt = (float)dl/30.0; // convert to time by dividing by c (30cm/ns) 
  return dt; 
}// end calcDeltaT

float skimmer::calcAvgT(const float t, const int n){ 
  float avgT = -1000;
  if (n!=0) avgT = (float)t/(float)n;
  return avgT;
}// end calcAvgT

float skimmer::smearVal(const float res){
  TRandom3 * r = new TRandom3(0);
  float val = r->Gaus(1,res);
  delete r;
  return val;
}// end smearVal

