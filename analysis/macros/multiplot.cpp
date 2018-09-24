#include "multiplot.hh"

multiplot::multiplot(TString cTau, TString inpath)
{

  ctau = cTau;
  path = inpath;
  odir = Form("%s/Compare/",path.Data()); 
  ofile = TFile::Open(Form("%scompare_ct%s.root",odir.Data(),ctau.Data()),"RECREATE");
  fin  = TFile::Open(Form("%s/XXqqqq_m50_ct%s/plots_dispjets.root",path.Data(),ctau.Data()),"READ");
  if (fin==NULL) std::cout << "File not found: " << fin << std::endl;

}// end multiplot constructor

multiplot::~multiplot()
{

  for (unsigned int i = 0; i < nh; i++){
    delete canv[i];
    delete leg[i];
    for (unsigned int j = 0; j < nsub[i]; j++){
      delete h[i][j];
    }
  }
  delete ofile;

}// end multiplot destructor

void multiplot::go()
{

  // setup histos to read
  multiplot::setuphistos();
  // setup cuts
  multiplot::setupcuts();
 
  float tmp_cut, tmp_int;
  float tmp_min, tmp_max, tmp_eff;
 
  h.resize(nh); // histograms
  e.resize(nh); // efficiencies
  for (unsigned int i = 0; i < nh; i++){ // loop over type of histos
    h[i].resize(nsub[i]);
    e[i].resize(nsub[i]);
    for (unsigned int j = 0; j < nsub[i]; j++){ // loop over individual histos

      // pick up the histograms
      h[i][j] = (TH1F*)fin->Get(Form("%s",hist[i][j].Data()));
      if (h[i][j]==NULL) std::cout << "NO HISTO" << std::endl;
      // get the efficiency
      e[i][j].resize(ncuts[i]);
      tmp_int = h[i][j]->Integral();
      tmp_max = h[i][j]->GetNbinsX(); // max bin 

      for (unsigned int k = 0; k < ncuts[i]; k++){ // loop over cuts 
        tmp_min = h[i][j]->FindBin(cuts[i][k]);
        tmp_cut = h[i][j]->Integral(tmp_min,tmp_max);
        e[i][j][k] = tmp_cut/tmp_int;
      }// end loop over cuts

      // rescale by integral  
      if (tmp_int > 0) h[i][j]->Scale(1/tmp_int); 

    }// end loop over individual histos 
  }// end loop over type of histos

  gStyle->SetOptStat(0);
  multiplot::prepcanvas();
  multiplot::drawplots();
  multiplot::saveplots();

  for (unsigned int i = 0; i < nh; i++){
    multiplot::setupeffplots(i);
  }

}// end go

void multiplot::setupeffplots(unsigned int i)
{

  // i = histogram type index
  // j = subhistogram index
  // k = cut index

  // setup canvas
  TCanvas * c = new TCanvas(Form("eff_ct%s_%s",ctau.Data(),hist[i][0].Data()),"",600,600);
  c->cd();

  // setup legend
  TLegend * l = new TLegend(0.46,0.66,0.86,0.86,NULL,"brNDC");
  l->SetLineColor(0);
  l->SetLineStyle(0);
  l->SetFillColor(0);
  l->SetFillStyle(0);

  // setup histo
  unsigned int nbin = cuts[i].size()-1;
  float cutarr[cuts[i].size()]; 
  std::copy(cuts[i].begin(), cuts[i].end(), cutarr);
  TH1FVec heff;
  heff.resize(nsub[i]);
  for (unsigned int j = 0; j < nsub[i]; j++){
    heff[j] = new TH1F(Form("heff_%i",j),"",nbin,cutarr);
    // fill for each cut
    for (unsigned int k = 0; k < ncuts[i]; k++){
      heff[j]->Fill(cuts[i][k],e[i][j][k]);
    }
  }

  // setup labels
  heff[0]->GetXaxis()->SetTitle("Cut value");
  heff[0]->GetYaxis()->SetTitle("Efficiency");

  // draw histo
  for (unsigned int j = 0; j < nsub[i]; j++){
    heff[j]->SetLineColor(colors[hist[i][j]]);
    if (j==0) heff[j]->Draw("HIST");
    else heff[j]->Draw("HIST SAME");
    l->AddEntry(heff[j],hist[i][j],"l"); 
  }
  l->Draw("SAME");

  // save
  c->SaveAs(Form("%sEff_ct%s_%s.png",odir.Data(),ctau.Data(),hist[i][0].Data()));
  c->SaveAs(Form("%sEff_ct%s_%s.pdf",odir.Data(),ctau.Data(),hist[i][0].Data()));
  ofile->cd();
  c->Write();
 
  // delete
  for (unsigned int j = 0; j < nsub[i]; j++) delete heff[j];
  delete l;
  delete c;
}

void multiplot::prepcanvas()
{

  canv.resize(nh);
  leg.resize(nh);
  for (unsigned int i = 0; i < nh; i++){
    canv[i] = new TCanvas(Form("%s_canv",hist[i][0].Data()),"",600,800);
    canv[i]->cd();

    leg[i] = new TLegend(0.5,0.7,0.9,0.9,NULL,"brNDC");
    leg[i]->SetLineColor(1);
    leg[i]->SetLineStyle(1);
    leg[i]->SetFillColor(0);
    leg[i]->SetFillStyle(0);

  }

}// end prepcanvas

void multiplot::drawplots()
{

  for (unsigned int i = 0; i < nh; i++){
    canv[i]->cd();
    for (unsigned int j = 0; j < nsub[i]; j++){
      h[i][j]->SetLineColor(colors[hist[i][j]]);
      if (j==0) h[i][j]->Draw("HIST");
      else      h[i][j]->Draw("HIST SAME");
      leg[i]->AddEntry(h[i][j],hist[i][j],"l"); 
    }
    leg[i]->Draw("SAME"); 
  } 

}// end drawplots

void multiplot::saveplots()
{

  for (unsigned int i = 0; i < nh; i++){
    canv[i]->cd();
    canv[i]->SetLogy(1);
    canv[i]->SaveAs(Form("%s%s_ct%s_log.png",odir.Data(),hist[i][0].Data(),ctau.Data()));
    canv[i]->SaveAs(Form("%s%s_ct%s_log.pdf",odir.Data(),hist[i][0].Data(),ctau.Data()));
    canv[i]->SetLogy(0);
    canv[i]->SaveAs(Form("%s%s_ct%s_lin.png",odir.Data(),hist[i][0].Data(),ctau.Data()));
    canv[i]->SaveAs(Form("%s%s_ct%s_lin.pdf",odir.Data(),hist[i][0].Data(),ctau.Data()));
    ofile->cd();
    canv[i]->Write(Form("%s",hist[i][0].Data())); 
  } 

}// end saveplots

void multiplot::setupcuts()
{

  cuts.resize(nh);

  for (unsigned int i = 0; i < 21; i++){
    cuts[0].push_back(i*0.1); // cuts[0] are for jet t plots
    cuts[1].push_back(i*0.1); // cuts[1] are for unmatch t plots
  }

  for (unsigned int i = 0; i < nh; i++){
    ncuts.push_back(cuts[i].size());
  }

}// end setupcuts

void multiplot::setuphistos()
{

  nh = 2;
  hist.resize(nh);
  // hist[0] are jet t plots
  hist[0].push_back("jet_t");
  hist[0].push_back("jet_t_smear30");
  hist[0].push_back("jet_t_smear50");
  hist[0].push_back("jet_t_smear70");
  hist[0].push_back("jet_t_smear180");
  hist[0].push_back("jet_t_smear500");
  // hist[1] are unmatched t plots
  hist[1].push_back("unmatch_t");
  hist[1].push_back("unmatch_t_smear30");
  hist[1].push_back("unmatch_t_smear50");
  hist[1].push_back("unmatch_t_smear70");
  hist[1].push_back("unmatch_t_smear180");
  hist[1].push_back("unmatch_t_smear500");

  for (unsigned int i = 0; i < nh; i++){
    nsub.push_back(hist[i].size());
  }

  colors["jet_t"] 		= kBlack;
  colors["jet_t_smear30"]	= kMagenta;
  colors["jet_t_smear50"]	= kTeal;
  colors["jet_t_smear70"]	= kGreen;
  colors["jet_t_smear180"]	= kRed;
  colors["jet_t_smear500"]	= kOrange;
  colors["unmatch_t"] 		= kBlack;
  colors["unmatch_t_smear30"]	= kMagenta;
  colors["unmatch_t_smear50"]	= kTeal;
  colors["unmatch_t_smear70"]	= kGreen;
  colors["unmatch_t_smear180"]	= kRed;
  colors["unmatch_t_smear500"]	= kOrange;

}// end setuphistos

