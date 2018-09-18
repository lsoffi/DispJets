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
 
  float tmp_int;
 
  h.resize(nh);
  for (unsigned int i = 0; i < nh; i++){
    h[i].resize(nsub[i]);
    for (unsigned int j = 0; j < nsub[i]; j++){
      // pick up the histograms
      h[i][j] = (TH1F*)fin->Get(Form("%s",hist[i][j].Data()));
      if (h[i][j]==NULL) std::cout << "NO HISTO" << std::endl;
      // rescale by integral  
      tmp_int = h[i][j]->Integral();
      if (tmp_int > 0) h[i][j]->Scale(1/tmp_int); 
    } 
  } 

  gStyle->SetOptStat(0);
  multiplot::prepcanvas();
  multiplot::drawplots();
  multiplot::saveplots();

}// end go

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
  // hist[1] are unmatched t plots
  hist[1].push_back("unmatch_t");
  hist[1].push_back("unmatch_t_smear30");
  hist[1].push_back("unmatch_t_smear50");
  hist[1].push_back("unmatch_t_smear70");
  hist[1].push_back("unmatch_t_smear180");

  for (unsigned int i = 0; i < nh; i++){
    nsub.push_back(hist[i].size());
  }

  colors["jet_t"] 		= kBlack;
  colors["jet_t_smear30"]	= kMagenta;
  colors["jet_t_smear50"]	= kTeal;
  colors["jet_t_smear70"]	= kGreen;
  colors["jet_t_smear180"]	= kRed;
  colors["unmatch_t"] 		= kBlack;
  colors["unmatch_t_smear30"]	= kMagenta;
  colors["unmatch_t_smear50"]	= kTeal;
  colors["unmatch_t_smear70"]	= kGreen;
  colors["unmatch_t_smear180"]	= kRed;

}// end setuphistos

