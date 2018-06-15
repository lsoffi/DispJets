import os
import sys
import ROOT
from array import array
import numpy as np
import json
from CMS_lumi import CMS_lumi
import plotting_interp as plot
ROOT.gROOT.SetBatch(ROOT.kTRUE)
plot.ModTDRStyle()

def run(model,m1,res,lumi):

  filepath = 'dispjets_jsons/'
  #m2  = [1, 3, 10, 30, 100, 300, 1000, 3000]
  m2 = [100]
  wgt = 1 

  # read in theoretical xsec
  if (model=="2HDM"): xsecfile = open('crosssectionZp2HDM.txt','r')
  if (model=="BARY"): xsecfile = open('crosssectionZpBaryonic.txt','r')
  if (model=="DISP"): xsecfile = open('crosssectionDispJets.txt','r') 
  mZ = []
  mA = []
  mstr = []
  xsecs = {} 
  xsec = []
  for line in xsecfile:
    line = line.rstrip('\n')
    line = line.split(' ')
    mZ.append(line[0])
    mA.append(line[1])
    massstr = line[0]+'_'+line[1]
    mstr.append(massstr)
    xsecs[massstr] = float(line[2])

  xsec = array( 'd' )
  mass = array( 'd' )
  exp_raw = array( 'd' ) 
  obs_raw = array( 'd' ) 
  up1_raw = array( 'd' )
  do1_raw = array( 'd' )
  up2_raw = array( 'd' )
  do2_raw = array( 'd' )

  # pick up jsons with limits
  for m in m2:
    newmassstr = str(m)+'_'+m1
    if newmassstr not in mstr: continue

    # make correct array of xsecs
    mass.append(m)
    xsec.append(xsecs[str(m)+'_'+m1])   
    if (model=="2HDM"): filename=filepath+'Zprime'+str(m)+'A'+m1+'.json'
    if (model=="BARY"): filename=filepath+'Zprime'+str(m)+'DM'+m1+'.json'
    if (model=="DISP"): filename=filepath+'XXto4Q_M'+m1+'_CT'+str(m)+'_r'+res+'_L'+lumi+'.json'

    if os.path.isfile(filename):
      with open(filename) as jsonfile:
        data = json.load(jsonfile)
        for key in data: 
          exp_raw.append(data[key][u'exp0'])
          obs_raw.append(data[key][u'obs'])
          up1_raw.append(data[key][u'exp+1'])
          do1_raw.append(data[key][u'exp-1'])
          up2_raw.append(data[key][u'exp+2'])
          do2_raw.append(data[key][u'exp-2'])
    else:  print 'File '+filename+' NOT found'

  exp = scaleBR(exp_raw,wgt) 
  obs = scaleBR(obs_raw,wgt)
  up1 = scaleBR(up1_raw,wgt)
  do1 = scaleBR(do1_raw,wgt)
  up2 = scaleBR(up2_raw,wgt)
  do2 = scaleBR(do2_raw,wgt)

  numpts = len(mass)          
  limitPlotExp = ROOT.TGraph(numpts,mass,exp)
  limitPlotObs = ROOT.TGraph(numpts,mass,obs)
  limitPlotXsc = ROOT.TGraph(numpts,mass,xsec)
  limitPlot1sig   = ROOT.TGraph(2*numpts)
  limitPlot2sig   = ROOT.TGraph(2*numpts)
  for i in range(0,numpts):
    limitPlot1sig.SetPoint(i,mass[i],up1[i])
    limitPlot1sig.SetPoint(numpts+i,mass[numpts-i-1],do1[numpts-i-1])
    limitPlot2sig.SetPoint(i,mass[i],up2[i])
    limitPlot2sig.SetPoint(numpts+i,mass[numpts-i-1],do2[numpts-i-1])

  ROOT.gStyle.SetOptStat(0)
  c = ROOT.TCanvas('','')
  c.SetLogx(1)
  c.SetLogy(1)
  c.SetGrid();
 
  limitPlotObs.SetTitle("")
  limitPlotObs.GetXaxis().SetTitle("LL particle c#tau [mm]")
  limitPlotObs.GetYaxis().SetTitle("#sigma [pb]")

  # styling
  limitPlotXsc.SetLineColor(4)
  limitPlotXsc.SetLineStyle(7)
  limitPlotXsc.SetLineWidth(2)

  limitPlotExp.SetLineStyle(7)
  limitPlotExp.SetLineWidth(2)
  limitPlotObs.SetLineWidth(2)

  limitPlot1sig.SetFillColor(5)
  limitPlot1sig.SetFillStyle(3013)

  limitPlot2sig.SetFillColor(8)
  limitPlot2sig.SetFillStyle(3013)

  limitPlotObs.SetMaximum(10)
  limitPlotObs.SetMinimum(0.1)

  limitPlotObs.Draw("ACP")
  limitPlot2sig.Draw("F SAME")
  limitPlot1sig.Draw("F SAME")
  limitPlotExp.Draw("CP SAME")
  limitPlotObs.Draw("CP SAME")
  limitPlotXsc.Draw("CP SAME")

  CMS_lumi(c,4,0)
  c.RedrawAxis() 
  c.Print("~/www/Plots/DispJets/GenLevelPlots/Limits/lim1D_m"+m1+"_r"+res+"_L"+lumi+".pdf")
  c.Print("~/www/Plots/DispJets/GenLevelPlots/Limits/lim1D_m"+m1+"_r"+res+"_L"+lumi+".png")
 
def scaleBR(val,wgt):
   scaledval = val
   return val
  
if __name__=="__main__":
  model = sys.argv[1]
  m1    = sys.argv[2]
  res   = sys.argv[3]
  lumi  = sys.argv[4]
  run(model,m1,res,lumi)
