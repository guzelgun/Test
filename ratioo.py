from ROOT import *
from utils import *
from histlib import *
import sys

gStyle.SetOptStat(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetLegendBorderSize(0)

try: fname = sys.argv[1]
except:
    fname = '/home/cakir/Work/Programs/CerenAnalysisCP/NtuplerAnalyzer/newVar4/SingleS_NoPU_483_300_his.root'


f  = TFile(fname)

keys = f.GetListOfKeys()

h  = f.Get("hAW1")
h2 = f.Get("hAW2")

overflow(h)
overflow(h2)
    
h3 = createRatio(h, h2)

h.SetLineWidth(2)
h.SetLineColor(kGreen +2)
h. SetLineStyle(1)

h2.SetLineWidth(2)
h2.SetLineColor(kRed)
h2. SetLineStyle(1)


c, pad1, pad2 = createCanvasPads()
c.cd()
leg = TLegend(0.60, 0.60, 0.87, 0.89)
leg.SetHeader("Signal point (cau, mass)")

leg.AddEntry(h, "AWH1 Distribution")
leg.AddEntry(h2, "AWH2 Distribution")
leg.SetFillStyle(0)
pad1.cd()
gPad.SetLogy()

h.GetXaxis().SetTitle("[GeV]")
    #h.GetXaxis().SetLabelSize(0.03)
    #h.GetYaxis().SetRangeUser(0.0001, (max(h.GetMaximum(),h2.GetMaximum())))
h.GetYaxis().SetTitle("Events norm to unity")
h.GetYaxis().SetTitleOffset(1)
h.SetTitle("AWH1 & AWH2 Distributions in SUSY Signal: 483_300")

h.Draw('hist')
h2.Draw("same")
leg.Draw()

h.GetYaxis().SetTitleSize(20)
h.GetYaxis().SetTitleFont(43)
h.GetYaxis().SetTitleOffset(1.55)
h.GetYaxis().SetLabelFont(43)
h.GetYaxis().SetLabelSize(15)
h.GetYaxis().SetLabelSize(15)
axis = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis.SetLabelFont(43)
axis.SetLabelSize(15)
axis.Draw()
pad2.cd()
h3.Draw("ep")
line = TLine(0, 1, 40, 1);
line.Draw()
pause()



