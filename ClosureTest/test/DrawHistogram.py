#!usr/bin/python
import ROOT
from ROOT import TFile, TH1F, THStack, TCanvas, TMath, gROOT, gPad, TLegend
#from setTDRStyle import setTDRStyle
from array import array

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--var",action="store",type="string",dest="var",default='ptHat')
parser.add_option("--xmin",action="store",type="float",dest="xmin",default=1)
parser.add_option("--xmax",action="store",type="float",dest="xmax",default=1)
parser.add_option("--xtitle",action="store",type="string",dest="xtitle",default='')
parser.add_option("--rebin",action="store",type="int",dest="rebin",default=1)
parser.add_option("--logy",action="store_true",default=False,dest="logy")
parser.add_option("--outputDir",action="store",type="string",default="./",dest="outputDir")

(options, args) = parser.parse_args()

var = options.var
xmin = options.xmin
xmax = options.xmax
xtitle = options.xtitle
rebin = options.rebin
logy = options.logy
outputDir = options.outputDir

# giulia -style
#gROOT.Reset()
#setTDRStyle()
#gROOT.ForceStyle()
#gROOT.SetStyle('tdrStyle')

fileNames = ['QCD_Pt-300to470','QCD_Pt-470to600','QCD_Pt-600to800', 'QCD_Pt-800to1000', 'QCD_Pt-1000to1400', 'QCD_Pt-1400to1800', 'QCD_Pt-1800to2400', 'QCD_Pt-2400to3200', 'QCD_Pt-3200', 'RSGravToJJ_kMpl01_M-1000', 'RSGravToJJ_kMpl01_M-5000', 'RSGravToJJ_kMpl01_M-8000']
xsections = [7475 ,587., 167, 28.25, 8.195, 0.7346, 0.102, 0.00644, 0.000163, 1, 1, 1]
colorF    = [ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9,ROOT.kBlue-9, ROOT.kWhite, ROOT.kWhite, ROOT.kWhite]
colorL    = [ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack,ROOT.kBlack, ROOT.kRed, ROOT.kGreen, ROOT.kBlue]
hist_noCuts      = []
hist_allCuts      = []
hist_allOtherCuts      = []
LUMI      = 1000.

## ---- CERN -------
PATH = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_0_pre9_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/data/output/'

#---- open the files --------------------
i_f = 0
for f in fileNames:
  inf = TFile.Open(PATH+'rootFile_'+f+'.root')
  print inf.GetName()
  
  Nev = inf.Get('cutHisto_noCuts_________________ptHat').GetEntries()
  wt = 1.0
  #if i_f < 3:
  wt = LUMI*xsections[i_f]/Nev

  h_noCuts = inf.Get('cutHisto_noCuts_________________'+var)
  h_noCuts.Scale(wt)
  h_noCuts.Rebin(rebin)
  h_noCuts.SetDirectory(0)
  h_noCuts.SetFillColor(colorF[i_f])
  h_noCuts.SetLineColor(colorL[i_f])
  h_noCuts.SetMarkerColor(colorL[i_f])
  hist_noCuts.append(h_noCuts)

  h_allOtherCuts = inf.Get('cutHisto_allOtherCuts___________'+var)
  h_allOtherCuts.Scale(wt)
  h_allOtherCuts.Rebin(rebin)
  h_allOtherCuts.SetDirectory(0)
  h_allOtherCuts.SetFillColor(colorF[i_f])
  h_allOtherCuts.SetLineColor(colorL[i_f])
  h_allOtherCuts.SetMarkerColor(colorL[i_f])
  hist_allOtherCuts.append(h_allOtherCuts)


  h_allCuts = inf.Get('cutHisto_allCuts________________'+var)
  h_allCuts.Scale(wt)
  h_allCuts.Rebin(rebin)
  h_allCuts.SetDirectory(0)
  h_allCuts.SetFillColor(colorF[i_f])
  h_allCuts.SetLineColor(colorL[i_f])
  h_allCuts.SetMarkerColor(colorL[i_f])
  hist_allCuts.append(h_allCuts)
   
  i_f += 1


#kFactor = NDAT/NQCD
kFactor = 1.3
print ("kFactor = %f" % kFactor)

                                      
for i in range(0,9) :
  hist_allCuts[i].Scale(kFactor)
  hist_allOtherCuts[i].Scale(kFactor)
  hist_noCuts[i].Scale(kFactor)


NQCD_allCuts = hist_allCuts[0].Integral()+hist_allCuts[1].Integral()+hist_allCuts[2].Integral()+hist_allCuts[3].Integral()+hist_allCuts[4].Integral()+hist_allCuts[5].Integral()+hist_allCuts[6].Integral()+hist_allCuts[7].Integral()+hist_allCuts[8].Integral()
NQCD_allOtherCuts = hist_allOtherCuts[0].Integral()+hist_allOtherCuts[1].Integral()+hist_allOtherCuts[2].Integral()+hist_allOtherCuts[3].Integral()+hist_allOtherCuts[4].Integral()+hist_allOtherCuts[5].Integral()+hist_allOtherCuts[6].Integral()+hist_allOtherCuts[7].Integral()+hist_allOtherCuts[8].Integral()
NQCD_noCuts = hist_noCuts[0].Integral()+hist_noCuts[1].Integral()+hist_noCuts[2].Integral()+hist_noCuts[3].Integral()+hist_noCuts[4].Integral()+hist_noCuts[5].Integral()+hist_noCuts[6].Integral()+hist_noCuts[7].Integral()+hist_noCuts[8].Integral()
     
hist_allCutsQCD = hist_allCuts[0].Clone('hist_allCutsQCD')
hist_noCutsQCD = hist_noCuts[0].Clone('hist_noCutsQCD')
hist_allOtherCutsQCD = hist_allOtherCuts[0].Clone('hist_allOtherCutsQCD')


for i in range(1,9):
  hist_allCutsQCD.Add(hist_allCuts[i])
  hist_allOtherCutsQCD.Add(hist_allOtherCuts[i])
  hist_noCutsQCD.Add(hist_noCuts[i])

hsQCD_allCuts = THStack('QCD_allCuts','QCD_allCuts')
hsQCD_allOtherCuts = THStack('QCD_allOtherCuts','QCD_allOtherCuts')
hsQCD_noCuts = THStack('QCD_noCuts','QCD_noCuts')

for i in range(0,9) :
  hsQCD_allCuts.Add(hist_allCuts[i])
  hsQCD_allOtherCuts.Add(hist_allOtherCuts[i])
  hsQCD_noCuts.Add(hist_noCuts[i])


hist_allCutsSig_1000 = hist_allCuts[9].Clone('hist_allCutsSig_1000')
hist_noCutsSig_1000 = hist_noCuts[9].Clone('hist_noCutsSig_1000')
hist_allOtherCutsSig_1000 = hist_allOtherCuts[9].Clone('hist_allOtherCutsSig_1000')
hist_allCutsSig_5000 = hist_allCuts[10].Clone('hist_allCutsSig_5000')
hist_noCutsSig_5000 = hist_noCuts[10].Clone('hist_noCutsSig_5000')
hist_allOtherCutsSig_5000 = hist_allOtherCuts[10].Clone('hist_allOtherCutsSig_5000')
hist_allCutsSig_8000 = hist_allCuts[11].Clone('hist_allCutsSig_8000')
hist_noCutsSig_8000 = hist_noCuts[11].Clone('hist_noCutsSig_8000')
hist_allOtherCutsSig_8000 = hist_allOtherCuts[11].Clone('hist_allOtherCutsSig_8000')


#signal scale as bkg if not mjj
#if (var!="Dijet_MassW") :# && var!="Dijet_MassAK8" && var!="Dijet_MassAK4") :
hist_noCutsSig_1000.Scale(NQCD_noCuts/hist_noCutsSig_1000.Integral())
hist_noCutsSig_5000.Scale(NQCD_noCuts/hist_noCutsSig_5000.Integral())
hist_noCutsSig_8000.Scale(NQCD_noCuts/hist_noCutsSig_8000.Integral())
hist_allOtherCutsSig_1000.Scale(NQCD_allOtherCuts/hist_allOtherCutsSig_1000.Integral())
hist_allOtherCutsSig_5000.Scale(NQCD_allOtherCuts/hist_allOtherCutsSig_5000.Integral())
hist_allOtherCutsSig_8000.Scale(NQCD_allOtherCuts/hist_allOtherCutsSig_8000.Integral())
hist_allCutsSig_1000.Scale(NQCD_allCuts/hist_allCutsSig_1000.Integral())
hist_allCutsSig_5000.Scale(NQCD_allCuts/hist_allCutsSig_5000.Integral())
hist_allCutsSig_8000.Scale(NQCD_allCuts/hist_allCutsSig_8000.Integral())

NSIG_allCuts_1000 = hist_allCutsSig_1000.Integral()
NSIG_allCuts_5000 = hist_allCutsSig_5000.Integral()
NSIG_allCuts_8000 = hist_allCutsSig_8000.Integral()

NSIG_allOtherCuts_1000 = hist_allOtherCutsSig_1000.Integral()
NSIG_allOtherCuts_5000 = hist_allOtherCutsSig_5000.Integral()
NSIG_allOtherCuts_8000 = hist_allOtherCutsSig_8000.Integral()

NSIG_noCuts_1000 = hist_noCutsSig_1000.Integral()
NSIG_noCuts_5000 = hist_noCutsSig_5000.Integral()
NSIG_noCuts_8000 = hist_noCutsSig_8000.Integral()



print ("---- After scaling signal to bkg (if not plotting mjj) -----")
print ("bkg integral no cuts = %f" % NQCD_noCuts)
print ("sig integral 1000 no cuts = %f" % NSIG_noCuts_1000)
print ("sig integral 5000 no cuts = %f" % NSIG_noCuts_5000)
print ("sig integral 8000 no cuts = %f" % NSIG_noCuts_8000)
print ("bkg integral allOther cuts = %f" % NQCD_allOtherCuts)
print ("sig integral 1000 allOther cuts = %f" % NSIG_allOtherCuts_1000)
print ("sig integral 5000 allOther cuts = %f" % NSIG_allOtherCuts_5000)
print ("sig integral 8000 allOther cuts = %f" % NSIG_allOtherCuts_8000)
print ("bkg integral all cuts = %f" % NQCD_allCuts)
print ("sig integral 1000 all cuts = %f" % NSIG_allCuts_1000)
print ("sig integral 5000 all cuts = %f" % NSIG_allCuts_5000)
print ("sig integral 8000 all cuts = %f" % NSIG_allCuts_8000)




#----- Drawing  and save on file -----------------------

outFile = TFile(outputDir+"histo_signal_bkg_"+var+".root", "recreate")
outFile.cd()
hist_allCutsQCD.Write()
hist_allOtherCutsQCD.Write()
hist_noCutsQCD.Write()
hist_allCutsSig_1000.Write()
hist_noCutsSig_1000.Write()
hist_allOtherCutsSig_1000.Write()
hist_allCutsSig_5000.Write()
hist_noCutsSig_5000.Write()
hist_allOtherCutsSig_5000.Write()
hist_allCutsSig_8000.Write()
hist_noCutsSig_8000.Write()
hist_allOtherCutsSig_8000.Write()


can_allCuts = TCanvas('can_allCuts_'+var,'can_allCuts_'+var,900,600)
can_allCuts_sig = TCanvas('can_allCuts_sig_'+var,'can_allCuts_sig_'+var,900,600)
can_allOtherCuts = TCanvas('can_allOtherCuts_'+var,'can_allOtherCuts_'+var,900,600)
can_allOtherCuts_sig = TCanvas('can_allOtherCuts_sig_'+var,'can_allOtherCuts_sig_'+var,900,600)
can_noCuts = TCanvas('can_noCuts_'+var,'can_noCuts_'+var,900,600)
can_noCuts_sig = TCanvas('can_noCuts_sig_'+var,'can_noCuts_sig_'+var,900,600)


leg = TLegend(0.6, 0.7, 0.85, 0.85)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.AddEntry(hist_allCutsQCD, "QCD", "l")
leg.AddEntry(hist_allCutsSig_1000, "RS Graviton M=1000", "l")
leg.AddEntry(hist_allCutsSig_5000, "RS Graviton M=5000", "l")
leg.AddEntry(hist_allCutsSig_8000, "RS Graviton M=8000", "l")


can_allCuts.cd()

if logy:
  gPad.SetLogy()

#hist_allCutsQCD.Reset()
hist_allCutsQCD.GetXaxis().SetRangeUser(xmin,xmax)
hist_allCutsQCD.GetXaxis().SetTitle(xtitle)
maximumBin = array('f',  [hist_allCutsQCD.GetBinContent(hist_allCutsQCD.GetMaximumBin()), hist_allCutsSig_1000.GetBinContent(hist_allCutsSig_1000.GetMaximumBin()), hist_allCutsSig_5000.GetBinContent(hist_allCutsSig_5000.GetMaximumBin()), hist_allCutsSig_8000.GetBinContent(hist_allCutsSig_8000.GetMaximumBin())])
max = TMath.MaxElement(4, maximumBin)
hist_allCutsQCD.SetMaximum(1.2*max)
hist_allCutsQCD.SetMinimum(0.00001)
hist_allCutsQCD.Draw("hist")
hist_allCutsSig_1000.Draw("hist same")
hist_allCutsSig_5000.Draw("hist same")
hist_allCutsSig_8000.Draw("hist same")
leg.Draw()
gPad.RedrawAxis()

can_allCuts.Write()   
can_allCuts.SaveAs(outputDir+var+'_allCuts.png')
can_allCuts.SaveAs(outputDir+var+'_allCuts.svg')
can_allCuts.Close()

can_allOtherCuts.cd()   

if logy:
  gPad.SetLogy()

#hist_allOtherCutsQCD.Reset()
hist_allOtherCutsQCD.GetXaxis().SetRangeUser(xmin,xmax)
hist_allOtherCutsQCD.GetXaxis().SetTitle(xtitle)
maximumBin = array('f',  [hist_allOtherCutsQCD.GetBinContent(hist_allOtherCutsQCD.GetMaximumBin()), hist_allOtherCutsSig_1000.GetBinContent(hist_allOtherCutsSig_1000.GetMaximumBin()), hist_allOtherCutsSig_5000.GetBinContent(hist_allOtherCutsSig_5000.GetMaximumBin()), hist_allOtherCutsSig_8000.GetBinContent(hist_allOtherCutsSig_8000.GetMaximumBin())])
max = TMath.MaxElement(4, maximumBin)
hist_allOtherCutsQCD.SetMaximum(1.2*max)

hist_allOtherCutsQCD.SetMinimum(0.00001)
hist_allOtherCutsQCD.Draw("hist")
hist_allOtherCutsSig_1000.Draw("hist same")
hist_allOtherCutsSig_5000.Draw("hist same")
hist_allOtherCutsSig_8000.Draw("hist same")
leg.Draw()
gPad.RedrawAxis()

can_allOtherCuts.Write()   
can_allOtherCuts.SaveAs(outputDir+var+'_allOtherCuts.png')
can_allOtherCuts.SaveAs(outputDir+var+'_allOtherCuts.svg')
can_allOtherCuts.Close()



can_noCuts.cd()   

if logy:
  gPad.SetLogy()

#hist_noCutsQCD.Reset()
hist_noCutsQCD.GetXaxis().SetRangeUser(xmin,xmax)
hist_noCutsQCD.GetXaxis().SetTitle(xtitle)

maximumBin = array('f',  [hist_noCutsQCD.GetBinContent(hist_noCutsQCD.GetMaximumBin()), hist_noCutsSig_1000.GetBinContent(hist_noCutsSig_1000.GetMaximumBin()), hist_noCutsSig_5000.GetBinContent(hist_noCutsSig_5000.GetMaximumBin()), hist_noCutsSig_8000.GetBinContent(hist_noCutsSig_8000.GetMaximumBin())])
max = TMath.MaxElement(4, maximumBin)
hist_noCutsQCD.SetMaximum(1.2*max)

hist_noCutsQCD.SetMinimum(0.00001)
hist_noCutsQCD.Draw("hist")
hist_noCutsSig_1000.Draw("hist same")
hist_noCutsSig_5000.Draw("hist same")
hist_noCutsSig_8000.Draw("hist same")
leg.Draw()

gPad.RedrawAxis()

can_noCuts.Write()   
can_noCuts.SaveAs(outputDir+var+'_noCuts.png')
can_noCuts.SaveAs(outputDir+var+'_noCuts.svg')
can_noCuts.Close()

##### signal only

leg_sig = TLegend(0.6, 0.7, 0.85, 0.85)
leg_sig.SetLineColor(0)
leg_sig.SetFillColor(0)
leg_sig.AddEntry(hist_allCutsSig_1000, "M=1000", "l")
leg_sig.AddEntry(hist_allCutsSig_5000, "M=5000", "l")
leg_sig.AddEntry(hist_allCutsSig_8000, "M=8000", "l")

can_allCuts_sig.cd()
hist_allCutsSig_1000.GetXaxis().SetRangeUser(xmin,xmax)
hist_allCutsSig_1000.GetXaxis().SetTitle(xtitle)
hist_allCutsSig_1000.GetYaxis().SetTitle("arb. units")
max_1000_5000 = TMath.Max(hist_allCutsSig_1000.GetBinContent(hist_allCutsSig_1000.GetMaximumBin()), hist_allCutsSig_5000.GetBinContent(hist_allCutsSig_5000.GetMaximumBin())) 
hist_allCutsSig_1000.SetMaximum(1.2*TMath.Max(hist_allCutsSig_8000.GetBinContent(hist_allCutsSig_8000.GetMaximumBin()), max_1000_5000))
hist_allCutsSig_1000.DrawNormalized("hist")
hist_allCutsSig_5000.DrawNormalized("hist same")
hist_allCutsSig_8000.DrawNormalized("hist same")
leg_sig.Draw()
can_allCuts_sig.Write()
can_allCuts_sig.SaveAs(outputDir+var+'_signal_allCuts.png')
can_allCuts_sig.SaveAs(outputDir+var+'_signal_allCuts.svg')
can_allCuts_sig.Close()


can_allOtherCuts_sig.cd()
hist_allOtherCutsSig_1000.GetXaxis().SetRangeUser(xmin,xmax)
hist_allOtherCutsSig_1000.GetXaxis().SetTitle(xtitle)
hist_allOtherCutsSig_1000.GetYaxis().SetTitle("arb. units")
max_1000_5000 = TMath.Max(hist_allOtherCutsSig_1000.GetBinContent(hist_allOtherCutsSig_1000.GetMaximumBin()), hist_allOtherCutsSig_5000.GetBinContent(hist_allOtherCutsSig_5000.GetMaximumBin())) 
hist_allOtherCutsSig_1000.SetMaximum(1.2*TMath.Max(hist_allOtherCutsSig_8000.GetBinContent(hist_allOtherCutsSig_8000.GetMaximumBin()), max_1000_5000))
hist_allOtherCutsSig_1000.DrawNormalized("hist")
hist_allOtherCutsSig_5000.DrawNormalized("hist same")
hist_allOtherCutsSig_8000.DrawNormalized("hist same")
leg_sig.Draw()
can_allOtherCuts_sig.Write()
can_allOtherCuts_sig.SaveAs(outputDir+var+'_signal_allOtherCuts.png')
can_allOtherCuts_sig.SaveAs(outputDir+var+'_signal_allOtherCuts.svg')
can_allOtherCuts_sig.Close()

can_noCuts_sig.cd()
hist_noCutsSig_1000.GetXaxis().SetRangeUser(xmin,xmax)
hist_noCutsSig_1000.GetXaxis().SetTitle(xtitle)
hist_noCutsSig_1000.GetYaxis().SetTitle("arb. units")
max_1000_5000 = TMath.Max(hist_noCutsSig_1000.GetBinContent(hist_noCutsSig_1000.GetMaximumBin()), hist_noCutsSig_5000.GetBinContent(hist_noCutsSig_5000.GetMaximumBin())) 
hist_noCutsSig_1000.SetMaximum(1.2*TMath.Max(hist_noCutsSig_8000.GetBinContent(hist_noCutsSig_8000.GetMaximumBin()), max_1000_5000))
hist_noCutsSig_1000.DrawNormalized("hist")
hist_noCutsSig_5000.DrawNormalized("hist same")
hist_noCutsSig_8000.DrawNormalized("hist same")
leg_sig.Draw()
can_noCuts_sig.Write()
can_noCuts_sig.SaveAs(outputDir+var+'_signal_noCuts.png')
can_noCuts_sig.SaveAs(outputDir+var+'_signal_noCuts.svg')
can_noCuts_sig.Close()

print ("hist_allCuts[9].GetEntries() = %d" % hist_allCuts[9].GetEntries())

outFile.Close()

#----- keep the GUI alive ------------
if __name__ == '__main__':
  rep = ''
  while not rep in ['q','Q']:
    rep = raw_input('enter "q" to quit: ')
    if 1 < len(rep):
      rep = rep[0]
