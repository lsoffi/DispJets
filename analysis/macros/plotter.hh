#ifndef __plotter__
#define __plotter__

#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
#include <vector>

typedef std::map< TString, TH1F*   > TH1map;
typedef std::map< TString, TH2F*   > TH2map;
typedef TH1map::iterator             TH1iter;
typedef TH2map::iterator             TH2iter;
typedef std::map< TString, TString > Samples;

class plotter
{
public:
  plotter(TString cTau, TString infile, TString outfile);
  ~plotter();
  void go();
  void setuptreebranches();
  void histos( TH1map & , TH2map & );
  void save1Dplots( const TH1map & );
  void save2Dplots( const TH2map & );
  void delete1Dplots( TH1map & );
  void delete2Dplots( TH2map & );
  TH1F * MakeTH1FPlot( const TString, const TString, const int, const double, const double, const TString, const TString );
  TH2F * MakeTH2FPlot( const TString, const TString, const int, const double, const double, const int, const double, const double, const TString, const TString );
  float calcDeltaT(const float lx, const float bx, const float la, const float ba, const float lo, const float bo);
  float smearVal(const float res);

private:
  Samples	infiles;
  int    	nsamples;
  TFile *       fin;
  TTree *       t;
  TString       odir;
  TFile *	fout;
  TString	ctau;

  // Declaration of leaf types
  Int_t           sample;
  Int_t           run;
  Int_t           lumi;
  Long64_t        event;
  Float_t         weight;
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
  vector<float>   *genpar_theta;
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
  vector<float>   *genpar_loline;
  vector<float>   *genpar_laline;
  vector<float>   *genpar_beta;
  vector<float>   *genpar_q;
  vector<int>     *genpar_match_q1lo;
  vector<int>     *genpar_match_q2lo;
  vector<int>     *genpar_match_q3lo;
  vector<int>     *genpar_match_q4lo;
  Int_t           nmothers;
  vector<int>     *mom_id;
  vector<int>     *mom_stat;
  vector<float>   *mom_e;
  vector<float>   *mom_m;
  vector<float>   *mom_p;
  vector<float>   *mom_pt;
  vector<float>   *mom_eta;
  vector<float>   *mom_phi;
  vector<float>   *mom_theta;
  vector<float>   *mom_vx;
  vector<float>   *mom_vy;
  vector<float>   *mom_vz;
  vector<float>   *mom_beta;
  vector<float>   *mom_gama;
  vector<float>   *mom_Lxy;
  vector<float>   *mom_Lz;
  vector<float>   *mom_Lxyz;
  vector<float>   *mom_ctau;
  vector<int>     *mom_dupl;
  vector<int>     *jet_nconst;
  vector<float>   *jet_avg_t;
  vector<float>   *jet_smear_30_t;
  vector<float>   *jet_smear_50_t;
  vector<float>   *jet_smear_70_t;
  vector<float>   *jet_smear_180_t;
  vector<float>   *jet_pt;
  vector<float>   *jet_alpha_PV;
  vector<float>   *jet_theta_2D;

  // List of branches
  TBranch        *b_sample;   //!
  TBranch        *b_run;   //!
  TBranch        *b_lumi;   //!
  TBranch        *b_event;   //!
  TBranch        *b_weight;   //!
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
  TBranch        *b_genpar_theta;   //!
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
  TBranch        *b_genpar_loline;   //!
  TBranch        *b_genpar_laline;   //!
  TBranch        *b_genpar_beta;   //!
  TBranch        *b_genpar_q;   //!
  TBranch        *b_genpar_match_q1lo;   //!
  TBranch        *b_genpar_match_q2lo;   //!
  TBranch        *b_genpar_match_q3lo;   //!
  TBranch        *b_genpar_match_q4lo;   //!
  TBranch        *b_nmothers;   //!
  TBranch        *b_mom_id;   //!
  TBranch        *b_mom_stat;   //!
  TBranch        *b_mom_e;   //!
  TBranch        *b_mom_m;   //!
  TBranch        *b_mom_p;   //!
  TBranch        *b_mom_pt;   //!
  TBranch        *b_mom_eta;   //!
  TBranch        *b_mom_phi;   //!
  TBranch        *b_mom_theta;   //!
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
  TBranch        *b_jet_nconst;   //!
  TBranch        *b_jet_avg_t;   //!
  TBranch        *b_jet_smear_30_t;   //!
  TBranch        *b_jet_smear_50_t;   //!
  TBranch        *b_jet_smear_70_t;   //!
  TBranch        *b_jet_smear_180_t;   //!
  TBranch        *b_jet_pt;   //!
  TBranch        *b_jet_alpha_PV;   //!
  TBranch        *b_jet_theta_2D;   //!


};

#endif
