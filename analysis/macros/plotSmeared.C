#include <iostream>

void plot( TString, TFile*, TString );

void plotSmeared(){

  TString path   = "~/www/Plots/DispJets/GenLevelPlots/NewSamples/";
  TString infile = path+"plots_dispjets.root";
  TFile * fin    = TFile::Open(infile);
  plot(path,fin,"jet_t");
  plot(path,fin,"unmatch_t"); 

}

void plot( TString path, TFile * fin, TString h ){
  
  gStyle->SetOptStat(0);

  TH1F *h0 = (TH1F*)fin->Get(Form("%s",h.Data()));
  TH1F *h1 = (TH1F*)fin->Get(Form("%s_smear30",h.Data()));
  TH1F *h2 = (TH1F*)fin->Get(Form("%s_smear180",h.Data()));

  h0->SetLineColor(kBlue);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kGreen);
  
  TCanvas *c = new TCanvas("c","c",600,800);
  c->cd();
  h0->Draw("HIST");
  h1->Draw("SAME HIST");
  //h2->Draw("SAME HIST");

  TLegend *l = new TLegend(0.5,0.7,0.9,0.8,NULL,"brNDC");
  l->SetBorderSize(0);
  l->AddEntry(h0,h,"L");
  l->AddEntry(h1,h+" smeared 30ns","L");
  //l->AddEntry(h2,h+" smeared 180ns","L");
  l->Draw("SAME");

  c->SetLogy(0);
  c->SaveAs(Form("%s%s_wSmear.png",path.Data(),h.Data()));
  c->SaveAs(Form("%s%s_wSmear.pdf",path.Data(),h.Data()));
  c->SetLogy(1);
  c->SaveAs(Form("%s%s_wSmear_log.png",path.Data(),h.Data()));
  c->SaveAs(Form("%s%s_wSmear_log.pdf",path.Data(),h.Data()));
   

}
