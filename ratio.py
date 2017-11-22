from ROOT import *
from utils import *
from histlib import *
import sys

gStyle.SetOptStat(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetLegendBorderSize(0)

try: fname = sys.argv[1]
except:
    fname = ''

try: fname2 = sys.argv[2]
except:
    fname2 = ''

f  = TFile(fname)
f2 = TFile(fname2)

keys = f.GetListOfKeys()

histlist = []

for key in keys:
    #c.SetLogy()
    hist = key.GetName()
    print hist
    histlist.append(hist[1:])
    h  = f.Get(hist)
    h2 = f2.Get(hist)
    overflow(h)
    overflow(h2)
    
    h3 = createRatio(h, h2)
#Normalise
    n = 1
    if h.Integral() > 0 and h2.Integral() > 0:
        s1 = n/(h.Integral())
        h.Scale(s1)
        s2 = n/(h2.Integral())
        h2.Scale(s2)

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

    leg.AddEntry(h, "Old Signal (1000 cm)")
    leg.AddEntry(h2, "New Signal (1000 cm)")
    leg.SetFillStyle(0)
    pad1.cd()
    gPad.SetLogy()

    h.GetXaxis().SetTitle("p_{T} of the 2^{nd} Hardest Jet")
    #h.GetXaxis().SetLabelSize(0.03)
    #h.GetYaxis().SetRangeUser(0.0001, (max(h.GetMaximum(),h2.GetMaximum())))
    h.GetYaxis().SetTitle("Events norm to unity")
    h.GetYaxis().SetTitleOffset(1)
    h.SetTitle(histlib[hist[1:]])

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
    pause()

print 'histlib= {}'
for hist in histlist:
    print 'histlib["' + hist + '"] = "' + hist +'"'



