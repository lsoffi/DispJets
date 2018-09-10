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
  vector<float>	*genpar_loline;		genpar_loline = 0;	TBranch *b_genpar_loline;
  vector<float>	*genpar_beta;		genpar_beta = 0;	TBranch	*b_genpar_beta; 
  vector<float>	*mom_beta;		mom_beta = 0;		TBranch *b_mom_beta;
  vector<float> *mom_Lxyz;		mom_Lxyz = 0;		TBranch *b_mom_Lxyz;
  vector<float> *mom_theta;		mom_theta = 0;		TBranch *b_mom_theta;
  vector<float> *genpar_theta;		genpar_theta = 0;	TBranch *b_genpar_theta;
  vector<int>	*genpar_match_q0;	genpar_match_q0 = 0;	TBranch *b_genpar_match_q0;
  vector<int>	*genpar_match_q1;	genpar_match_q1 = 0;	TBranch *b_genpar_match_q1;
  vector<int>	*genpar_match_q2;	genpar_match_q2 = 0;	TBranch *b_genpar_match_q2;
  vector<int>	*genpar_match_q3;	genpar_match_q3 = 0;	TBranch *b_genpar_match_q3;
  vector<int>	*genpar_match_q4;	genpar_match_q4 = 0;	TBranch *b_genpar_match_q4;
  vector<int>	*genpar_match_q1lo;	genpar_match_q1lo = 0;	TBranch *b_genpar_match_q1lo;
  vector<int>	*genpar_match_q2lo;	genpar_match_q2lo = 0;	TBranch *b_genpar_match_q2lo;
  vector<int>	*genpar_match_q3lo;	genpar_match_q3lo = 0;	TBranch *b_genpar_match_q3lo;
  vector<int>	*genpar_match_q4lo;	genpar_match_q4lo = 0;	TBranch *b_genpar_match_q4lo;
  
  intree->SetBranchAddress("ngenpart", 		&ngenpart, 		&b_ngenpart);
  intree->SetBranchAddress("genpar_stat", 	&genpar_stat, 		&b_genpar_stat);
  intree->SetBranchAddress("genpar_eta", 	&genpar_eta, 		&b_genpar_eta);
  intree->SetBranchAddress("genpar_pt", 	&genpar_pt,		&b_genpar_pt);
  intree->SetBranchAddress("genpar_la", 	&genpar_la,		&b_genpar_la);
  intree->SetBranchAddress("genpar_loline", 	&genpar_loline, 	&b_genpar_loline);
  intree->SetBranchAddress("genpar_beta", 	&genpar_beta, 		&b_genpar_beta);
  intree->SetBranchAddress("mom_beta", 		&mom_beta, 		&b_mom_beta);
  intree->SetBranchAddress("mom_Lxyz", 		&mom_Lxyz, 		&b_mom_Lxyz);
  intree->SetBranchAddress("mom_theta",		&mom_theta,		&b_mom_theta);
  intree->SetBranchAddress("genpar_theta",	&genpar_theta,		&b_genpar_theta);
  intree->SetBranchAddress("genpar_match_q0",	&genpar_match_q0,	&b_genpar_match_q0); 
  intree->SetBranchAddress("genpar_match_q1",	&genpar_match_q1,	&b_genpar_match_q1); 
  intree->SetBranchAddress("genpar_match_q2",	&genpar_match_q2,	&b_genpar_match_q2); 
  intree->SetBranchAddress("genpar_match_q3",	&genpar_match_q3,	&b_genpar_match_q3); 
  intree->SetBranchAddress("genpar_match_q4",	&genpar_match_q4,	&b_genpar_match_q4); 
  intree->SetBranchAddress("genpar_match_q1lo",	&genpar_match_q1lo,	&b_genpar_match_q1lo); 
  intree->SetBranchAddress("genpar_match_q2lo",	&genpar_match_q2lo,	&b_genpar_match_q2lo); 
  intree->SetBranchAddress("genpar_match_q3lo",	&genpar_match_q3lo,	&b_genpar_match_q3lo); 
  intree->SetBranchAddress("genpar_match_q4lo",	&genpar_match_q4lo,	&b_genpar_match_q4lo); 

  // variables to add to the tree
  vector<float> genpar_t;		outtree->Branch("genpar_t",		&genpar_t);
  vector<int>   jet_nconst;		outtree->Branch("jet_nconst",		&jet_nconst);
  vector<float> jet_avg_t;		outtree->Branch("jet_avg_t",		&jet_avg_t);
  vector<float> jet_smear_30_t;		outtree->Branch("jet_smear_30_t",	&jet_smear_30_t);
  vector<float> jet_smear_50_t;		outtree->Branch("jet_smear_50_t",	&jet_smear_50_t);
  vector<float> jet_smear_70_t;		outtree->Branch("jet_smear_70_t",	&jet_smear_70_t);
  vector<float> jet_smear_180_t;	outtree->Branch("jet_smear_180_t",	&jet_smear_180_t);
  vector<float> jet_pt;			outtree->Branch("jet_pt",		&jet_pt);
  vector<float> jet_alpha_PV;		outtree->Branch("jet_alpha_PV",		&jet_alpha_PV);
  vector<float> jet_theta_2D;		outtree->Branch("jet_theta_2D",		&jet_theta_2D);

  // loop over events
  for (unsigned int entry = 0; entry < intree->GetEntries(); entry++){
    intree->GetEntry(entry);

    // mother info
    float jet1_mom_beta  = (*mom_beta)[0];
    float jet1_mom_lxyz  = (*mom_Lxyz)[0];
    float jet2_mom_beta  = (*mom_beta)[1];
    float jet2_mom_lxyz  = (*mom_Lxyz)[1];
    float jet3_mom_beta  = (*mom_beta)[2];
    float jet3_mom_lxyz  = (*mom_Lxyz)[2];
    float jet4_mom_beta  = (*mom_beta)[3];
    float jet4_mom_lxyz  = (*mom_Lxyz)[3];
    float jet1_mom_theta = (*mom_theta)[0];
    float jet2_mom_theta = (*mom_theta)[1];
    float jet3_mom_theta = (*mom_theta)[2];
    float jet4_mom_theta = (*mom_theta)[3];

    vector<float> tmp_genpar_t;
    vector<int>   tmp_jet_nconst;
    vector<float> tmp_jet_avg_t;	
    vector<float> tmp_jet_smear_30_t;	
    vector<float> tmp_jet_smear_50_t;	
    vector<float> tmp_jet_smear_70_t;	
    vector<float> tmp_jet_smear_180_t;
    vector<float> tmp_jet_pt;	
    vector<float> tmp_jet_alpha_PV;
    vector<float> tmp_jet_theta_2D;

    int   jet1_nconst         = 0;
    float jet1_pt             = 0;
    float jet1_pt_all         = 0;
    float jet1_const_dt       = 0; 
    float jet1_time_raw       = 0;
    float jet1_time_smear_30  = 0; 
    float jet1_time_smear_50  = 0; 
    float jet1_time_smear_70  = 0; 
    float jet1_time_smear_180 = 0;
    int   jet2_nconst         = 0; 
    float jet2_pt             = 0;
    float jet2_pt_all         = 0;
    float jet2_const_dt       = 0; 
    float jet2_time_raw       = 0;
    float jet2_time_smear_30  = 0; 
    float jet2_time_smear_50  = 0; 
    float jet2_time_smear_70  = 0; 
    float jet2_time_smear_180 = 0;
    int   jet3_nconst         = 0; 
    float jet3_pt             = 0;
    float jet3_pt_all         = 0;
    float jet3_const_dt       = 0; 
    float jet3_time_raw       = 0;
    float jet3_time_smear_30  = 0; 
    float jet3_time_smear_50  = 0; 
    float jet3_time_smear_70  = 0; 
    float jet3_time_smear_180 = 0;
    int   jet4_nconst         = 0; 
    float jet4_pt             = 0;
    float jet4_pt_all         = 0;
    float jet4_const_dt       = 0; 
    float jet4_time_raw       = 0;
    float jet4_time_smear_30  = 0; 
    float jet4_time_smear_50  = 0; 
    float jet4_time_smear_70  = 0; 
    float jet4_time_smear_180 = 0;

    float smear_value_30      = 0; 
    float smear_value_50      = 0; 
    float smear_value_70      = 0; 
    float smear_value_180     = 0;
    float jet_const_beta      = 0; 
    float jet_const_lxyz      = 0;
    float jet_orig_beta       = 0; 
    float jet_orig_lxyz       = 0;  
 
    float tmp_const_pt_jet1    = 0;
    float tmp_const_pt_jet2    = 0;
    float tmp_const_pt_jet3    = 0;
    float tmp_const_pt_jet4    = 0;
    float tmp_const_theta_jet1 = 0;
    float tmp_const_theta_jet2 = 0;
    float tmp_const_theta_jet3 = 0;
    float tmp_const_theta_jet4 = 0;


    // loop over gen particles
    for (int gp = 0; gp < ngenpart; gp++){

      if ((*genpar_stat)[gp] != 1)            continue; // keep only final state particles
      if (std::fabs((*genpar_eta)[gp]) > 1.5) continue; // keep only particles in the barrel
      if ((*genpar_pt)[gp] < 1.0)             continue; // keep only pT > 1 GeV tracks
  
      smear_value_30  = smearVal(0.03); // res: 30 ps -> 0.03 ns
      smear_value_50  = smearVal(0.05); // res: 50 ps -> 0.05 ns
      smear_value_70  = smearVal(0.07); // res: 70 ps -> 0.07 ns
      smear_value_180 = smearVal(0.18); // res: 180 ps -> 0.18 ns
 
      // constituent particle info 
      jet_const_beta = 1.0;
      jet_const_lxyz = (*genpar_la)[gp];
      // hypothetical non-LL info
      jet_orig_beta  = 1.0; // (*genpar_beta)[gp]
      jet_orig_lxyz  = (*genpar_loline)[gp];    
     

      // sum jet pt of loose matched particles
      if ( (*genpar_match_q1lo)[gp]==1 ) jet1_pt_all += (*genpar_pt)[gp];
      if ( (*genpar_match_q2lo)[gp]==1 ) jet2_pt_all += (*genpar_pt)[gp];
      if ( (*genpar_match_q3lo)[gp]==1 ) jet3_pt_all += (*genpar_pt)[gp];
      if ( (*genpar_match_q4lo)[gp]==1 ) jet4_pt_all += (*genpar_pt)[gp];
 
      //--- q1 jet
      if ( (*genpar_match_q1)[gp]==1 ){
        jet1_const_dt        = calcDeltaT(jet1_mom_lxyz,jet1_mom_beta,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta);
        jet1_nconst         += 1.0; 
        jet1_time_raw       += jet1_const_dt;
        jet1_time_smear_30  += jet1_const_dt*smear_value_30;
        jet1_time_smear_50  += jet1_const_dt*smear_value_50;
        jet1_time_smear_70  += jet1_const_dt*smear_value_70;
        jet1_time_smear_180 += jet1_const_dt*smear_value_180;
        jet1_pt             += (*genpar_pt)[gp]; 
        if ( (*genpar_pt)[gp] >= tmp_const_pt_jet1 ){ // get theta of highest momentum constituent
          tmp_const_pt_jet1    = (*genpar_pt)[gp];
          tmp_const_theta_jet1 = (*genpar_theta)[gp]; 
        }
      }// end match to q1 jet

      // --- q2 jet
      if ( (*genpar_match_q2)[gp]==1){
        jet2_const_dt        = calcDeltaT(jet2_mom_lxyz,jet2_mom_beta,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta);
        jet2_nconst         += 1.0;
        jet2_time_raw       += jet2_const_dt;
        jet2_time_smear_30  += jet2_const_dt*smear_value_30;
        jet2_time_smear_50  += jet2_const_dt*smear_value_50;
        jet2_time_smear_70  += jet2_const_dt*smear_value_70;
        jet2_time_smear_180 += jet2_const_dt*smear_value_180;
        jet2_pt             += (*genpar_pt)[gp]; 
        if ( (*genpar_pt)[gp] >= tmp_const_pt_jet2 ){ // get theta of highest momentum constituent
          tmp_const_pt_jet2    = (*genpar_pt)[gp];
          tmp_const_theta_jet2 = (*genpar_theta)[gp]; 
        }
      }// end match to q2 jet

      // --- q3 jet
      if ( (*genpar_match_q3)[gp]==1){
        jet3_const_dt        = calcDeltaT(jet3_mom_lxyz,jet3_mom_beta,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta);
        jet3_nconst         += 1.0;
        jet3_time_raw       += jet3_const_dt;
        jet3_time_smear_30  += jet3_const_dt*smear_value_30;
        jet3_time_smear_50  += jet3_const_dt*smear_value_50;
        jet3_time_smear_70  += jet3_const_dt*smear_value_70;
        jet3_time_smear_180 += jet3_const_dt*smear_value_180;
        jet3_pt             += (*genpar_pt)[gp]; 
        if ( (*genpar_pt)[gp] >= tmp_const_pt_jet3 ){ // get theta of highest momentum constituent
          tmp_const_pt_jet3    = (*genpar_pt)[gp];
          tmp_const_theta_jet3 = (*genpar_theta)[gp]; 
        }
      }// end match to q3 jet

      // --- q4 jet
      if ( (*genpar_match_q4)[gp]==1){
        jet4_const_dt        = calcDeltaT(jet4_mom_lxyz,jet4_mom_beta,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta);
        jet4_nconst         += 1.0;
        jet4_time_raw       += jet4_const_dt;
        jet4_time_smear_30  += jet4_const_dt*smear_value_30;
        jet4_time_smear_50  += jet4_const_dt*smear_value_50;
        jet4_time_smear_70  += jet4_const_dt*smear_value_70;
        jet4_time_smear_180 += jet4_const_dt*smear_value_180;
        jet4_pt             += (*genpar_pt)[gp]; 
        if ( (*genpar_pt)[gp] >= tmp_const_pt_jet4 ){ // get theta of highest momentum constituent
          tmp_const_pt_jet4    = (*genpar_pt)[gp];
          tmp_const_theta_jet4 = (*genpar_theta)[gp]; 
        }
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
    // alpha jet variable
    tmp_jet_alpha_PV.push_back((jet1_pt_all-jet1_pt)/jet1_pt_all);
    tmp_jet_alpha_PV.push_back((jet2_pt_all-jet2_pt)/jet2_pt_all);
    tmp_jet_alpha_PV.push_back((jet3_pt_all-jet3_pt)/jet3_pt_all);
    tmp_jet_alpha_PV.push_back((jet4_pt_all-jet4_pt)/jet4_pt_all);
    // theta 2D variable
    tmp_jet_theta_2D.push_back((jet1_mom_theta-tmp_const_theta_jet1));
    tmp_jet_theta_2D.push_back((jet2_mom_theta-tmp_const_theta_jet2));
    tmp_jet_theta_2D.push_back((jet3_mom_theta-tmp_const_theta_jet3));
    tmp_jet_theta_2D.push_back((jet4_mom_theta-tmp_const_theta_jet4));
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
    // smeared jet time (res 50 ps)
    tmp_jet_smear_50_t.push_back(calcAvgT(jet1_time_smear_50,jet1_nconst));
    tmp_jet_smear_50_t.push_back(calcAvgT(jet2_time_smear_50,jet2_nconst));
    tmp_jet_smear_50_t.push_back(calcAvgT(jet3_time_smear_50,jet3_nconst));
    tmp_jet_smear_50_t.push_back(calcAvgT(jet4_time_smear_50,jet4_nconst));
    // smeared jet time (res 70 ps)
    tmp_jet_smear_70_t.push_back(calcAvgT(jet1_time_smear_70,jet1_nconst));
    tmp_jet_smear_70_t.push_back(calcAvgT(jet2_time_smear_70,jet2_nconst));
    tmp_jet_smear_70_t.push_back(calcAvgT(jet3_time_smear_70,jet3_nconst));
    tmp_jet_smear_70_t.push_back(calcAvgT(jet4_time_smear_70,jet4_nconst));
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
    jet_smear_50_t.clear();
    jet_smear_70_t.clear();
    jet_smear_180_t.clear();  
    jet_alpha_PV.clear();
    jet_theta_2D.clear();

    for (int i = 0; i < 4; i++){
      jet_nconst.push_back(tmp_jet_nconst[i]);
      jet_pt.push_back(tmp_jet_pt[i]);
      jet_avg_t.push_back(tmp_jet_avg_t[i]);
      jet_smear_30_t.push_back(tmp_jet_smear_30_t[i]);
      jet_smear_50_t.push_back(tmp_jet_smear_50_t[i]);
      jet_smear_70_t.push_back(tmp_jet_smear_70_t[i]);
      jet_smear_180_t.push_back(tmp_jet_smear_180_t[i]);
      jet_alpha_PV.push_back(tmp_jet_alpha_PV[i]);
      jet_theta_2D.push_back(tmp_jet_theta_2D[i]);
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

