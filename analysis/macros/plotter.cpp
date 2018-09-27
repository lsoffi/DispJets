#include "plotter.hh"

plotter::plotter(TString cTau, TString infile, TString outpath)
{

  ctau = cTau;
  odir = outpath;
  fin  = TFile::Open(infile,"READ"); 
  if (fin == NULL) std::cout << "File not found: " << infile << std::endl;
  fout = TFile::Open(Form("%splots_dispjets.root",odir.Data()),"RECREATE");

}// end plotter constructor

plotter::~plotter()
{

  std::cout << "Finishing up" << std::endl;
  delete t;
  delete fin;
  delete fout;

}// end plotter destructor

void plotter::histos( TH1map & map , TH2map & map2){

  // 1d histos
  map["LL_cTau"] = MakeTH1FPlot("LL_cTau","",2000,0,2000,"LL particle c#tau [cm]","");
  //if (ctau=="0mm") 		map["LL_cTau"] = MakeTH1FPlot("LL_cTau","",50,0,5,"LL particle c#tau [cm]","");
  //if (ctau=="1mm") 		map["LL_cTau"] = MakeTH1FPlot("LL_cTau","",50,0,5,"LL particle c#tau [cm]","");
  //if (ctau=="10mm") 		map["LL_cTau"] = MakeTH1FPlot("LL_cTau","",50,0,5,"LL particle c#tau [cm]","");
  //if (ctau=="100mm") 		map["LL_cTau"] = MakeTH1FPlot("LL_cTau","",500,0,50,"LL particle c#tau [cm]","");
  //if (ctau=="1000mm") 		map["LL_cTau"] = MakeTH1FPlot("LL_cTau","",500,0,500,"LL particle c#tau [cm]","");
  //if (ctau=="10000mm") 		map["LL_cTau"] = MakeTH1FPlot("LL_cTau","",5000,0,5000,"LL particle c#tau [cm]","");
  map["LL_beta"]		= MakeTH1FPlot("LL_beta","",100,0,1,"LL particle #beta","");
  map["particle_dxy"]		= MakeTH1FPlot("particle_dxy","",100,0,100,"d_{xy} [cm]","");
  map["particle_dxyz"]		= MakeTH1FPlot("particle_dxyz","",100,0,100,"d_{xyz} [cm]","");
  map["jet_alpha"]		= MakeTH1FPlot("jet_alpha","",100,0,1,"Jet #alpha","");
  map["jet_theta2D"]		= MakeTH1FPlot("jet_theta2D","",70,-5,2,"Jet #theta_{2D}","");
  map["nconst"]			= MakeTH1FPlot("nconst","",50,0,50,"Num. jet constituents","");
  map["const_t"]		= MakeTH1FPlot("const_t","",150,-5,10,"Jet constituent time [ns]",""); 
  map["const_pt"]		= MakeTH1FPlot("const_pt","",100,0,100,"Jet constituent p_{T} [GeV]","");
  map["max_jet_t"]		= MakeTH1FPlot("max_jet_t","",150,-5,10,"Max jet time [ns]","");
  map["jet_t"]			= MakeTH1FPlot("jet_t","",150,-5,10,"Jet time [ns]",""); 
  map["jet_t_smear0"]		= MakeTH1FPlot("jet_t_smear0","",150,-5,10,"Jet time [ns]","");
  map["jet_t_smear30"]		= MakeTH1FPlot("jet_t_smear30","",150,-5,10,"Jet time [ns]","");
  map["jet_t_smear50"]		= MakeTH1FPlot("jet_t_smear50","",150,-5,10,"Jet time [ns]","");
  map["jet_t_smear70"]		= MakeTH1FPlot("jet_t_smear70","",150,-5,10,"Jet time [ns]","");
  map["jet_t_smear180"] 	= MakeTH1FPlot("jet_t_smear180","",150,-5,10,"Jet time [ns]","");
  map["jet_t_smear500"] 	= MakeTH1FPlot("jet_t_smear500","",150,-5,10,"Jet time [ns]","");
  map["jet_t_diff0"]		= MakeTH1FPlot("jet_t_diff0","",200,-10,10,"T_{smeared} - T_{gen} [ns]","");
  map["jet_t_diff30"]		= MakeTH1FPlot("jet_t_diff30","",200,-10,10,"T_{smeared} - T_{gen} [ns]","");
  map["jet_t_diff50"]		= MakeTH1FPlot("jet_t_diff50","",200,-10,10,"T_{smeared} - T_{gen} [ns]","");
  map["jet_t_diff70"]		= MakeTH1FPlot("jet_t_diff70","",200,-10,10,"T_{smeared} - T_{gen} [ns]","");
  map["jet_t_diff180"]		= MakeTH1FPlot("jet_t_diff180","",200,-10,10,"T_{smeared} - T_{gen} [ns]","");
  map["jet_t_diff500"]		= MakeTH1FPlot("jet_t_diff500","",200,-10,10,"T_{smeared} - T_{gen} [ns]","");
  map["unmatch_beta"]		= MakeTH1FPlot("unmatch_beta","",100,0,1,"Unmatched particle #beta","");
  map["unmatch_t"]		= MakeTH1FPlot("unmatch_t","",150,-5,10,"Unmatched particle time [ns]","");
  map["unmatch_t_smear30"]	= MakeTH1FPlot("unmatch_t_smear30","",150,-5,10,"Unmatched particle time [ns]","");
  map["unmatch_t_smear50"]	= MakeTH1FPlot("unmatch_t_smear50","",150,-5,10,"Unmatched particle time [ns]","");
  map["unmatch_t_smear70"]	= MakeTH1FPlot("unmatch_t_smear70","",150,-5,10,"Unmatched particle time [ns]","");
  map["unmatch_t_smear180"]	= MakeTH1FPlot("unmatch_t_smear180","",150,-5,10,"Unmatched particle time [ns]","");
  map["unmatch_t_smear500"]	= MakeTH1FPlot("unmatch_t_smear500","",150,-5,10,"Unmatched particle time [ns]","");
  map["genjet_match_alpha"]	= MakeTH1FPlot("genjet_match_alpha","",100,0,1,"Jet #alpha","");
  map["genjet_match_theta"]	= MakeTH1FPlot("genjet_match_theta","",70,-5,2,"Jet #theta_{2D}","");
  map["genjet_unmatch_alpha"]	= MakeTH1FPlot("genjet_unmatch_alpha","",100,0,1,"Jet #alpha","");
  map["genjet_unmatch_theta"]	= MakeTH1FPlot("genjet_unmatch_theta","",70,-5,2,"Jet #theta_{2D}","");

  //// 2d histos
  //map2["tX1_tX2"]		= MakeTH2FPlot("tX1_tX2","",150,-5,10,150,-5,10,"Time jets from X1 [ns]","Time jets from X2 [ns]");
  //map2["tnear_tfar"]		= MakeTH2FPlot("tnear_tfar","",150,-5,10,250,-5,20,"Time close jets [ns]","Time far jets [ns]");

}

void plotter::go(){

  // get input tree
  t = (TTree*)fin->Get("dispjets/tree");
  if (t == NULL) std::cout << "Tree not found" << std::endl;
  setuptreebranches();

  // setup histo map
  TH1map h1map;
  TH2map h2map;
  histos(h1map, h2map);

   //--- loop over events
   unsigned int nentries = t->GetEntries();
   for (unsigned int entry = 0; entry < nentries; entry++){
     t->GetEntry(entry);

     // draw mom information
     for (unsigned int gm = 0; gm < 4; gm++){
       if ((*mom_dupl)[gm]==1) continue; // only plot mom once
       h1map["LL_beta"]->Fill((*mom_beta)[gm]);          // mom beta
       h1map["LL_cTau"]->Fill((*mom_ctau)[gm]);          // mom ctau
     }

     // draw gen jet information
     for (unsigned int gj = 0; gj < 4; gj++){
       if ( std::fabs((*genjet_eta)[gj]) > 1.5 ) continue; // keep only jets in barrel

       if ( (*genjet_match)[gj]==1 ){ // jet matches a LL particle
         h1map["genjet_match_alpha"]->Fill((*genjet_alpha_PV)[gj]);
         h1map["genjet_match_theta"]->Fill((*genjet_theta_2D)[gj]);
       }
       else { // jet does not match a LL particle
         h1map["genjet_unmatch_alpha"]->Fill((*genjet_alpha_PV)[gj]);
         h1map["genjet_unmatch_theta"]->Fill((*genjet_theta_2D)[gj]);
       }
     }

     // draw constituent information
     for (int gp = 0; gp < ngenpart; gp++){

         if ( (*genpar_stat)[gp] != 1 ) continue;           // keep only final state particles
         if ( std::abs((*genpar_eta)[gp]) > 1.5 ) continue; // keep only particles in barrel
         if ( (*genpar_pt)[gp] < 1.0 ) continue;            // keep only tracks of pT > 1 GeV 

         // constituent particle info
         float jet_const_beta = 1.0;
         float jet_const_lxyz = (*genpar_la)[gp];
         // hypothetical non-LL info
         float jet_orig_beta  = 1.0;
         float jet_orig_lxyz  = (*genpar_loline)[gp];     

         h1map["particle_dxy"]->Fill((*genpar_Lxy)[gp]);
         h1map["particle_dxyz"]->Fill((*genpar_Lxyz)[gp]);

         if ( (*genpar_match_q1)[gp]==1 || (*genpar_match_q2)[gp]==1 || (*genpar_match_q3)[gp]==1 || (*genpar_match_q4)[gp]==1 ){
             h1map["const_pt"]->Fill((*genpar_pt)[gp]);
             float jet1_const_t = calcDeltaT((*mom_Lxyz)[0],(*mom_beta)[0],jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta); 
             float jet2_const_t = calcDeltaT((*mom_Lxyz)[1],(*mom_beta)[1],jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta); 
             float jet3_const_t = calcDeltaT((*mom_Lxyz)[2],(*mom_beta)[2],jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta); 
             float jet4_const_t = calcDeltaT((*mom_Lxyz)[3],(*mom_beta)[3],jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta); 
             if ( (*genpar_match_q1)[gp]==1 ) h1map["const_t"]->Fill(jet1_const_t);
             if ( (*genpar_match_q2)[gp]==1 ) h1map["const_t"]->Fill(jet2_const_t);
             if ( (*genpar_match_q3)[gp]==1 ) h1map["const_t"]->Fill(jet3_const_t);
             if ( (*genpar_match_q4)[gp]==1 ) h1map["const_t"]->Fill(jet4_const_t);
         }// end jet constituent conditional

         if ( (*genpar_match_q0)[gp]==1 ){
             h1map["unmatch_beta"]->Fill((*genpar_beta)[gp]);
             float unmatch_t = calcDeltaT(0.0,1.0,jet_const_lxyz,jet_const_beta,jet_orig_lxyz,jet_orig_beta); 
             h1map["unmatch_t"]->Fill(unmatch_t);
             h1map["unmatch_t_smear30"]->Fill(unmatch_t*smearVal(0.03));
             h1map["unmatch_t_smear50"]->Fill(unmatch_t*smearVal(0.05));
             h1map["unmatch_t_smear70"]->Fill(unmatch_t*smearVal(0.07));
             h1map["unmatch_t_smear180"]->Fill(unmatch_t*smearVal(0.18));
             h1map["unmatch_t_smear500"]->Fill(unmatch_t*smearVal(0.50));
         }// end unmatched conditional

     }// end loop over genparticles 

     float max_jet_t = 0;
     // draw jet information (calculated in the skim)
     for (int i = 0; i < 4; i++){
         h1map["nconst"]->Fill((*jet_nconst)[i]);              // number of constituents in each jet
         h1map["jet_alpha"]->Fill((*jet_alpha_PV)[i]);         // alpha PV of each jet
         h1map["jet_theta2D"]->Fill((*jet_theta_2D)[i]);       // theta 2D of each jet
         h1map["jet_t"]->Fill((*jet_avg_t)[i]);                // jet (averaged constituent) timing
         h1map["jet_t_smear0"]->Fill((*jet_smear_0_t)[i]);     // smeared 0ps histos
         h1map["jet_t_smear30"]->Fill((*jet_smear_30_t)[i]);   // smeared 30ps histos
         h1map["jet_t_smear50"]->Fill((*jet_smear_50_t)[i]);   // smeared 50ps histos
         h1map["jet_t_smear70"]->Fill((*jet_smear_70_t)[i]);   // smeared 70ps histos
         h1map["jet_t_smear180"]->Fill((*jet_smear_180_t)[i]); // smeared 180ps histos
         h1map["jet_t_smear500"]->Fill((*jet_smear_500_t)[i]); // smeared 500ps histos
         h1map["jet_t_diff0"]->Fill(((*jet_smear_0_t)[i]-(*jet_avg_t)[i])); // diff smeared and jet time
         h1map["jet_t_diff30"]->Fill(((*jet_smear_30_t)[i]-(*jet_avg_t)[i])); // diff smeared and jet time
         h1map["jet_t_diff50"]->Fill(((*jet_smear_50_t)[i]-(*jet_avg_t)[i])); // diff smeared and jet time
         h1map["jet_t_diff70"]->Fill(((*jet_smear_70_t)[i]-(*jet_avg_t)[i])); // diff smeared and jet time
         h1map["jet_t_diff180"]->Fill(((*jet_smear_180_t)[i]-(*jet_avg_t)[i])); // diff smeared and jet time
         h1map["jet_t_diff500"]->Fill(((*jet_smear_500_t)[i]-(*jet_avg_t)[i])); // diff smeared and jet time
         if ((*jet_avg_t)[i] > max_jet_t){ h1map["max_jet_t"]->Fill(max_jet_t); max_jet_t = (*jet_avg_t)[i]; } // max avg jet time
     }

     //float dt_jetsX1   = 0;
     //float dt_jetsX2   = 0;
     //float dt_jet_near = 0;
     //float dt_jet_far  = 0;
     //h2map["tX1_tX2"]->Fill(dt_jetsX1,dt_jetsX2);
     //h2map["tnear_tfar"]->Fill(dt_jet_near,dt_jet_far);

  }// end loop over events

  // fit ctau plot
  h1map["LL_cTau"]->Fit("expo");

  // fit time diff plots
  h1map["jet_t_diff0"]->Fit("gaus");
  h1map["jet_t_diff30"]->Fit("gaus");
  h1map["jet_t_diff50"]->Fit("gaus");
  h1map["jet_t_diff70"]->Fit("gaus");
  h1map["jet_t_diff180"]->Fit("gaus");
  h1map["jet_t_diff500"]->Fit("gaus");

  // save plots
  std::cout << "Saving plots" << std::endl;
  save1Dplots(h1map);
  save2Dplots(h2map);

  // delete histos
  delete1Dplots(h1map);
  delete2Dplots(h2map);

}// end go()


float plotter::calcDeltaT(const float lx, const float bx, const float la, const float ba, const float lo, const float bo){
  float dl = (float)lx/(float)bx + (float)la/(float)ba - (float)lo/(float)bo; // change in distance
  float dt = (float)dl/30.0; // convert to time by dividing by c (30cm/ns) 
  return dt; 
}// end calcDeltaT

float plotter::smearVal(const float res){
  TRandom3 * r = new TRandom3(0);
  float val = r->Gaus(1,res);
  delete r;
  return val;
}// end smearVal

TH1F * plotter::MakeTH1FPlot(const TString hname, const TString htitle, const int nbins, const double xlow, const double xhigh, const TString xtitle, const TString ytitle)
{
  TString ytitleNew;
  Float_t binwidth = (xhigh-xlow)/nbins;
  if (ytitle=="") ytitleNew = "Events"; //Form("Events / %2.1f GeV",binwidth);
  else ytitleNew = ytitle;
 
  TH1F * hist = new TH1F(hname.Data(),htitle.Data(),nbins,xlow,xhigh);
  hist->GetXaxis()->SetTitle(xtitle.Data());
  hist->GetYaxis()->SetTitle(ytitleNew.Data());
  hist->Sumw2();
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit(1);
  return hist;
}// end MakeTH1FPlot

TH2F * plotter::MakeTH2FPlot(const TString hname, const TString htitle, const Int_t xnbins, const Double_t xlow, const Double_t xhigh, const Int_t ynbins, const Double_t ylow, const Double_t yhigh, const TString xtitle, const TString ytitle)
{
  TH2F * hist = new TH2F(hname.Data(),htitle.Data(),xnbins,xlow,xhigh,ynbins,ylow,yhigh);
  hist->GetXaxis()->SetTitle(xtitle.Data());
  hist->GetYaxis()->SetTitle(ytitle.Data());
  return hist;
}// end MakeTH2FPlot

void plotter::save1Dplots(const TH1map & map)
{
  fout->cd();
  std::cout << "Saving" << std::endl; 

  for (const auto & h:map){ // loop over histos
    const auto & name = h.first;
    const auto & hist = h.second;

    // save to output file
    hist->Write(hist->GetName(),TObject::kWriteDelete);

    // draw to canvas
    TCanvas *c = new TCanvas("c","c",600,800);
    c->cd();
    hist->Draw("HIST");

    // specific instructions for certain histos
    if (name=="LL_cTau"){
      TF1* f = (TF1*)hist->GetFunction("expo");
      f->Draw("SAME");
    }

    c->SetLogy(0);
    c->SaveAs(Form("%s%s.png",odir.Data(),name.Data()));
    c->SaveAs(Form("%s%s.pdf",odir.Data(),name.Data()));
    c->SetLogy(1);
    c->SaveAs(Form("%s%s_log.png",odir.Data(),name.Data()));
    c->SaveAs(Form("%s%s_log.pdf",odir.Data(),name.Data()));

    delete c;
  }// end loop over histos

}// end save1Dplots

void plotter::save2Dplots(const TH2map & map)
{
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

void plotter::delete1Dplots(TH1map & map)
{
  for ( TH1iter iter = map.begin(); iter != map.end(); ++iter ){
    delete (iter->second);
  }
  map.clear(); 
}// end delete1Dplots

void plotter::delete2Dplots(TH2map & map)
{
  for ( TH2iter iter = map.begin(); iter != map.end(); ++iter ){
    delete (iter->second);
  }
  map.clear(); 
}// end delete2Dplots

void plotter::setuptreebranches()
{
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
   genjet_match = 0;
   genjet0_const_st = 0;
   genjet0_const_id = 0;
   genjet0_const_pt = 0;
   genjet0_const_pv = 0;
   genjet0_const_theta = 0;
   genjet0_const_mtheta = 0;
   genjet0_const_vx = 0;
   genjet0_const_vy = 0;
   genjet0_const_vz = 0;
   genjet1_const_st = 0;
   genjet1_const_id = 0;
   genjet1_const_pt = 0;
   genjet1_const_pv = 0;
   genjet1_const_theta = 0;
   genjet1_const_mtheta = 0;
   genjet1_const_vx = 0;
   genjet1_const_vy = 0;
   genjet1_const_vz = 0;
   genjet2_const_st = 0;
   genjet2_const_id = 0;
   genjet2_const_pt = 0;
   genjet2_const_pv = 0;
   genjet2_const_theta = 0;
   genjet2_const_mtheta = 0;
   genjet2_const_vx = 0;
   genjet2_const_vy = 0;
   genjet2_const_vz = 0;
   genjet3_const_st = 0;
   genjet3_const_id = 0;
   genjet3_const_pt = 0;
   genjet3_const_pv = 0;
   genjet3_const_theta = 0;
   genjet3_const_mtheta = 0;
   genjet3_const_vx = 0;
   genjet3_const_vy = 0;
   genjet3_const_vz = 0;
   genpar_id = 0;
   genpar_stat = 0;
   genpar_pt = 0;
   genpar_eta = 0;
   genpar_phi = 0;
   genpar_theta = 0;
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
   genpar_dv1 = 0;
   genpar_dv2 = 0;
   genpar_dv3 = 0;
   genpar_dv4 = 0;
   genpar_lo = 0;
   genpar_la = 0;
   genpar_loline = 0;
   genpar_laline = 0;
   genpar_beta = 0;
   genpar_q = 0;
   genpar_match_q1lo = 0;
   genpar_match_q2lo = 0;
   genpar_match_q3lo = 0;
   genpar_match_q4lo = 0;
   mom_id = 0;
   mom_stat = 0;
   mom_e = 0;
   mom_m = 0;
   mom_p = 0;
   mom_pt = 0;
   mom_eta = 0;
   mom_phi = 0;
   mom_theta = 0;
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
   genpar_t = 0;
   jet_nconst = 0;
   jet_avg_t = 0;
   jet_smear_0_t = 0;
   jet_smear_30_t = 0;
   jet_smear_50_t = 0;
   jet_smear_70_t = 0;
   jet_smear_180_t = 0;
   jet_smear_500_t = 0;
   jet_pt = 0;
   jet_alpha_PV = 0;
   jet_theta_2D = 0;
   genjet_alpha_PV = 0;
   genjet_theta_2D = 0;
  
   // set branches
   t->SetBranchAddress("sample", &sample, &b_sample);
   t->SetBranchAddress("run", &run, &b_run);
   t->SetBranchAddress("lumi", &lumi, &b_lumi);
   t->SetBranchAddress("event", &event, &b_event);
   t->SetBranchAddress("weight", &weight, &b_weight);
   t->SetBranchAddress("ngenjets", &ngenjets, &b_ngenjets);
   t->SetBranchAddress("genjet0_nconst", &genjet0_nconst, &b_genjet0_nconst);
   t->SetBranchAddress("genjet1_nconst", &genjet1_nconst, &b_genjet1_nconst);
   t->SetBranchAddress("genjet2_nconst", &genjet2_nconst, &b_genjet2_nconst);
   t->SetBranchAddress("genjet3_nconst", &genjet3_nconst, &b_genjet3_nconst);
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
   t->SetBranchAddress("genjet0_const_pv", &genjet0_const_pv, &b_genjet0_const_pv);
   t->SetBranchAddress("genjet0_const_theta", &genjet0_const_theta, &b_genjet0_const_theta);
   t->SetBranchAddress("genjet0_const_mtheta", &genjet0_const_mtheta, &b_genjet0_const_mtheta);
   t->SetBranchAddress("genjet0_const_vx", &genjet0_const_vx, &b_genjet0_const_vx);
   t->SetBranchAddress("genjet0_const_vy", &genjet0_const_vy, &b_genjet0_const_vy);
   t->SetBranchAddress("genjet0_const_vz", &genjet0_const_vz, &b_genjet0_const_vz);
   t->SetBranchAddress("genjet1_const_st", &genjet1_const_st, &b_genjet1_const_st);
   t->SetBranchAddress("genjet1_const_id", &genjet1_const_id, &b_genjet1_const_id);
   t->SetBranchAddress("genjet1_const_pt", &genjet1_const_pt, &b_genjet1_const_pt);
   t->SetBranchAddress("genjet1_const_pv", &genjet1_const_pv, &b_genjet1_const_pv);
   t->SetBranchAddress("genjet1_const_theta", &genjet1_const_theta, &b_genjet1_const_theta);
   t->SetBranchAddress("genjet1_const_mtheta", &genjet1_const_mtheta, &b_genjet1_const_mtheta);
   t->SetBranchAddress("genjet1_const_vx", &genjet1_const_vx, &b_genjet1_const_vx);
   t->SetBranchAddress("genjet1_const_vy", &genjet1_const_vy, &b_genjet1_const_vy);
   t->SetBranchAddress("genjet1_const_vz", &genjet1_const_vz, &b_genjet1_const_vz);
   t->SetBranchAddress("genjet2_const_st", &genjet2_const_st, &b_genjet2_const_st);
   t->SetBranchAddress("genjet2_const_id", &genjet2_const_id, &b_genjet2_const_id);
   t->SetBranchAddress("genjet2_const_pt", &genjet2_const_pt, &b_genjet2_const_pt);
   t->SetBranchAddress("genjet2_const_pv", &genjet2_const_pv, &b_genjet2_const_pv);
   t->SetBranchAddress("genjet2_const_theta", &genjet2_const_theta, &b_genjet2_const_theta);
   t->SetBranchAddress("genjet2_const_mtheta", &genjet2_const_mtheta, &b_genjet2_const_mtheta);
   t->SetBranchAddress("genjet2_const_vx", &genjet2_const_vx, &b_genjet2_const_vx);
   t->SetBranchAddress("genjet2_const_vy", &genjet2_const_vy, &b_genjet2_const_vy);
   t->SetBranchAddress("genjet2_const_vz", &genjet2_const_vz, &b_genjet2_const_vz);
   t->SetBranchAddress("genjet3_const_st", &genjet3_const_st, &b_genjet3_const_st);
   t->SetBranchAddress("genjet3_const_id", &genjet3_const_id, &b_genjet3_const_id);
   t->SetBranchAddress("genjet3_const_pt", &genjet3_const_pt, &b_genjet3_const_pt);
   t->SetBranchAddress("genjet3_const_pv", &genjet3_const_pv, &b_genjet3_const_pv);
   t->SetBranchAddress("genjet3_const_theta", &genjet3_const_theta, &b_genjet3_const_theta);
   t->SetBranchAddress("genjet3_const_mtheta", &genjet3_const_mtheta, &b_genjet3_const_mtheta);
   t->SetBranchAddress("genjet3_const_vx", &genjet3_const_vx, &b_genjet3_const_vx);
   t->SetBranchAddress("genjet3_const_vy", &genjet3_const_vy, &b_genjet3_const_vy);
   t->SetBranchAddress("genjet3_const_vz", &genjet3_const_vz, &b_genjet3_const_vz);
   t->SetBranchAddress("ngenpart", &ngenpart, &b_ngenpart);
   t->SetBranchAddress("genpar_id", &genpar_id, &b_genpar_id);
   t->SetBranchAddress("genpar_stat", &genpar_stat, &b_genpar_stat);
   t->SetBranchAddress("genpar_pt", &genpar_pt, &b_genpar_pt);
   t->SetBranchAddress("genpar_eta", &genpar_eta, &b_genpar_eta);
   t->SetBranchAddress("genpar_phi", &genpar_phi, &b_genpar_phi);
   t->SetBranchAddress("genpar_theta", &genpar_theta, &b_genpar_theta);
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
   t->SetBranchAddress("genpar_dv1", &genpar_dv1, &b_genpar_dv1);
   t->SetBranchAddress("genpar_dv2", &genpar_dv2, &b_genpar_dv2);
   t->SetBranchAddress("genpar_dv3", &genpar_dv3, &b_genpar_dv3);
   t->SetBranchAddress("genpar_dv4", &genpar_dv4, &b_genpar_dv4);
   t->SetBranchAddress("genpar_lo", &genpar_lo, &b_genpar_lo);
   t->SetBranchAddress("genpar_la", &genpar_la, &b_genpar_la);
   t->SetBranchAddress("genpar_loline", &genpar_loline, &b_genpar_loline);
   t->SetBranchAddress("genpar_laline", &genpar_laline, &b_genpar_laline);
   t->SetBranchAddress("genpar_beta", &genpar_beta, &b_genpar_beta);
   t->SetBranchAddress("genpar_q", &genpar_q, &b_genpar_q);
   t->SetBranchAddress("genpar_match_q1lo", &genpar_match_q1lo, &b_genpar_match_q1lo);
   t->SetBranchAddress("genpar_match_q2lo", &genpar_match_q2lo, &b_genpar_match_q2lo);
   t->SetBranchAddress("genpar_match_q3lo", &genpar_match_q3lo, &b_genpar_match_q3lo);
   t->SetBranchAddress("genpar_match_q4lo", &genpar_match_q4lo, &b_genpar_match_q4lo);
   t->SetBranchAddress("nmothers", &nmothers, &b_nmothers);
   t->SetBranchAddress("mom_id", &mom_id, &b_mom_id);
   t->SetBranchAddress("mom_stat", &mom_stat, &b_mom_stat);
   t->SetBranchAddress("mom_e", &mom_e, &b_mom_e);
   t->SetBranchAddress("mom_m", &mom_m, &b_mom_m);
   t->SetBranchAddress("mom_p", &mom_p, &b_mom_p);
   t->SetBranchAddress("mom_pt", &mom_pt, &b_mom_pt);
   t->SetBranchAddress("mom_eta", &mom_eta, &b_mom_eta);
   t->SetBranchAddress("mom_phi", &mom_phi, &b_mom_phi);
   t->SetBranchAddress("mom_theta", &mom_theta, &b_mom_theta);
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
   t->SetBranchAddress("genpar_t", &genpar_t, &b_genpar_t);
   t->SetBranchAddress("jet_nconst", &jet_nconst, &b_jet_nconst);
   t->SetBranchAddress("jet_avg_t", &jet_avg_t, &b_jet_avg_t);
   t->SetBranchAddress("jet_smear_0_t", &jet_smear_0_t, &b_jet_smear_0_t);
   t->SetBranchAddress("jet_smear_30_t", &jet_smear_30_t, &b_jet_smear_30_t);
   t->SetBranchAddress("jet_smear_50_t", &jet_smear_50_t, &b_jet_smear_50_t);
   t->SetBranchAddress("jet_smear_70_t", &jet_smear_70_t, &b_jet_smear_70_t);
   t->SetBranchAddress("jet_smear_180_t", &jet_smear_180_t, &b_jet_smear_180_t);
   t->SetBranchAddress("jet_smear_500_t", &jet_smear_500_t, &b_jet_smear_500_t);
   t->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   t->SetBranchAddress("jet_alpha_PV", &jet_alpha_PV, &b_jet_alpha_PV);
   t->SetBranchAddress("jet_theta_2D", &jet_theta_2D, &b_jet_theta_2D);
   t->SetBranchAddress("genjet_alpha_PV", &genjet_alpha_PV, &b_genjet_alpha_PV);
   t->SetBranchAddress("genjet_theta_2D", &genjet_theta_2D, &b_genjet_theta_2D);

}
