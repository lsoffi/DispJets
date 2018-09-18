#include "overlay.hh"

overlay::overlay(std::vector< TString > cTau, TString inpath, Names files, Colors colormap)
{

  ctau = cTau;
  path = inpath;
  nfiles = ctau.size();
  for (unsigned int i = 0; i < nfiles; i++){
    infile[ctau[i]] = TFile::Open(files[ctau[i]],"READ");
    if (infile[ctau[i]]==NULL) std::cout << "File not found: " << files[ctau[i]] <<  std::endl;
  }
  odir = Form("%s/Overlay/",path.Data());
  ofile = TFile::Open(Form("%splots_overlay.root",odir.Data()),"RECREATE"); 
  colors = colormap;

} // end overlay constructor

overlay::~overlay()
{

  for (unsigned int i = 0; i < nh; i++){
    delete canv[i];
    delete leg[i];
    for (unsigned int j = 0; j < nfiles; j++){
      delete histos[i][j];
      delete ohistos[i][j];
    }
  }

  delete ofile;

}// end overlay destructor

void overlay::go()
{

  // setup names of plots to overlay for all samples
  overlay::setupplots();

  float tmp_int;
  float h_max = 0;
  float tmp_max = 0;
 
  // get the histos
  histos.resize(nh);
  ohistos.resize(nh);
  max_val.resize(nh);
  for (unsigned int i = 0; i < nh; i++){
    histos[i].resize(nfiles);
    ohistos[i].resize(nfiles);
    for (unsigned int j = 0; j < nfiles; j++){
      // pick up histos
      histos[i][j] = (TH1F*)infile[ctau[j]]->Get(Form("%s",hnames[i].Data()));
      if (histos[i][j]==NULL) std::cout << "NO HISTO" << std::endl;
      ohistos[i][j] = (TH1F*)histos[i][j]->Clone();
      // rescale by integral
      tmp_int = histos[i][j]->Integral();
      if (tmp_int > 0 ) ohistos[i][j]->Scale(1/tmp_int);
      // get maximum value of the histo
      h_max   = ohistos[i][j]->GetMaximum();
      tmp_max = (h_max >= tmp_max)? h_max : tmp_max;  
    }
    max_val[i] = tmp_max; 
  }

  gStyle->SetOptStat(0);
  overlay::prepcanvas();
  overlay::drawplots();
  overlay::saveplots();

}// end go

void overlay::drawplots()
{

  for (unsigned int i = 0; i < nh; i++){
    canv[i]->cd();
    for (unsigned int j = 0; j < nfiles; j++){
      ohistos[i][j]->SetLineColor(colors[ctau[j]]);
      if (j==0){
        ohistos[i][j]->SetMaximum(max_val[i]*10);
        ohistos[i][j]->Draw("AXIS");
      }
      else ohistos[i][j]->Draw("HIST SAME");
      leg[i]->AddEntry(ohistos[i][j],TString::Format("%s",ctau[j].Data()),"l");
    }
    leg[i]->Draw("SAME");
  } 

}// end drawplots

void overlay::saveplots()
{

  for (unsigned int i = 0; i < nh; i++){
    canv[i]->cd();
    canv[i]->SetLogy(1);
    canv[i]->SaveAs(Form("%s%s_log.png",odir.Data(),hnames[i].Data()));
    canv[i]->SaveAs(Form("%s%s_log.pdf",odir.Data(),hnames[i].Data()));
    canv[i]->SetLogy(0);
    canv[i]->SaveAs(Form("%s%s_lin.png",odir.Data(),hnames[i].Data()));
    canv[i]->SaveAs(Form("%s%s_lin.pdf",odir.Data(),hnames[i].Data()));
    ofile->cd();
    canv[i]->Write(Form("%s",hnames[i].Data())); 
  } 

}// end saveplots

void overlay::prepcanvas()
{

  // canvas
  canv.resize(nh);
  // legend
  leg.resize(nh);
  for (unsigned int i = 0; i < nh; i++){
    canv[i] = new TCanvas(Form("%s_canv",hnames[i].Data()),"",600,800);
    canv[i]->cd();

    leg[i] = new TLegend(0.5,0.7,0.9,0.9,NULL,"brNDC");
    leg[i]->SetLineColor(1);
    leg[i]->SetLineStyle(1);
    leg[i]->SetFillColor(0);
    leg[i]->SetFillStyle(0);

  } 

}// end prepcanvas

void overlay::setupplots()
{

  hnames.push_back("jet_alpha");
  hnames.push_back("jet_theta2D");
  hnames.push_back("jet_t");
  hnames.push_back("jet_t_smear30");
  nh = hnames.size();

}// end setupplots

