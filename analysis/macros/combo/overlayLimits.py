import sys
from ROOT import *
from optparse import OptionParser

def run(opts):

  # input files
  channels = [] 
  channels.append("new")                 # our new limits
  if opts.do_th: channels.append("thr")  # limits from theorists
  if opts.do_r2: channels.append("old")  # limits from runII analysis

  filepath = {}
  filepath["new"] = "dispjets_jsons/plot_m50_r30_L3000000.root"
  filepath["thr"] = "compare_files/data_lim_m50_theory.txt"
  filepath["old"] = "compare_files/data_lim_m50_L2600_inpb_dividedBy13p19.txt"

  # style for plots
  color = {}
  color["new"] = kRed
  color["thr"] = kGreen+1
  color["old"] = kBlue-4
  text = {}
  text["new"] = "w/ MTD (30ps) (3 ab-1)"
  text["thr"] = "arXiv:1805.05957v3  (3 ab-1)"
  text["old"] = "Run2 EXO-16-003 (2.6 fb-1)"

  # pickup graphs
  tgraph = {} 
  for ch in channels:
    if ch=="new":
      tgraph[ch] = TFile(filepath[ch]).Get("expected_curve") 
    else:
      tgraph[ch] = TGraph(filepath[ch])
 
  # change theory x from m -> mm
  rescaleaxis(tgraph["thr"],"x",1e3) 


  # setup drawing
  c = TCanvas("","",800,800)
  c.SetLogx(1)
  c.SetLogy(1)
  c.SetTicks()
  c.cd()
  for ch in channels:
    tgraph[ch].SetTitle("")
    tgraph[ch].GetXaxis().SetTitle("LL particle c#tau [mm]")
    tgraph[ch].GetYaxis().SetTitle("UL 95% CL on BR(H #rightarrow XX)")
    tgraph[ch].SetMaximum(2e3)
    tgraph[ch].SetLineColor(color[ch])
    tgraph[ch].SetLineStyle(1)
    tgraph[ch].SetLineWidth(2)
  tgraph["thr"].Draw("LA")
  tgraph["old"].Draw("L SAME") # limits at like 10^3
  tgraph["new"].Draw("L SAME")

  l = TLegend(0.5,0.6,0.85,0.85)
  l.SetFillStyle(0)
  l.SetBorderSize(0)
  for ch in channels: l.AddEntry(tgraph[ch],text[ch],"L")
  l.Draw("SAME")

  # save 
  suffix = ""
  if opts.suffix!="": suffix = "_"+opts.suffix 
  c.SaveAs(opts.outdir+"limits_compare"+suffix+".png") 
  c.SaveAs(opts.outdir+"limits_compare"+suffix+".pdf") 


def rescaleaxis(g,axis,scale):
    """This function rescales the x-axis on a TGraph."""
    N = g.GetN()
    if axis=="x": x = g.GetX()
    if axis=="y": x = g.GetY()
    for i in range(N):
        x[i] *= scale
    g.GetHistogram().Delete()
    g.SetHistogram(0)
    return

def init():
  # options
  parser = OptionParser("usage: %prog [options]")
  parser.add_option("-o","--outdir",action="store",type="string",
                    default="~/www/DispJets/Limits/",help="Output directory [default = %default]"), 
  parser.add_option("--suffix",action="store",type="string",
                    default="",help="String for output plot name [default = %default]"),
  parser.add_option("--theory",action="store_true",dest="do_th",
                    default=False,help="Plot limits from theorists [default = %default]"),
  parser.add_option("--orig",action="store_true",dest="do_r2",
                    default=False,help="Plot limits from RunII results [default = %default]"),
  (options, args) = parser.parse_args()
  
  # run
  run(options)

if __name__=="__main__":
  init() 
