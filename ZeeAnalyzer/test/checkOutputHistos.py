import ROOT


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

infile=ROOT.TFile.Open("outfileEleComparisonDY_noThr_15GeV.root","READ")

pRawMass=infile.Get("pInvMassRaw")
rRawMass=infile.Get("rInvMassRaw")
print 'number of entries %i' %pRawMass.GetEntries()
nbins=pRawMass.GetXaxis().GetNbins()
for ibin in range(0, nbins):
    pEntries=pRawMass.GetBinContent(ibin)
    rEntries=rRawMass.GetBinContent(ibin)
    print "bin %i, thr=0 number of entries=%i, thr=15 GeV number of entries=%i" %(ibin, pEntries, rEntries)
    if (pEntries-rEntries)!=0: print "difference: %i" %(pEntries-rEntries)

canvas=ROOT.TCanvas("c", "c", 600, 600)
canvas.Draw()
pad1=ROOT.TPad("pad1", "pad1", 0.,0.33,1,1)
pad2=ROOT.TPad("pad2", "pad2", 0.,0.,1,0.3)
pad1.Draw()
pad2.Draw()
pad1.cd()
pRawMass.Draw()
rRawMass.Draw("SAME")
rRawMass.SetLineColor(2)
pRawMass.GetXaxis().SetTitle("GeV")
leg=ROOT.TLegend(0.15,0.85, 0.35, 0.7)
leg.AddEntry(pRawMass, "thr=0GeV", "l")
leg.AddEntry(rRawMass, "thr=15GeV", "l")
leg.Draw()

pad2.cd()
#pRawMass.Sumw2()
#rRawMass.Sumw2()
ratio=pRawMass.Clone()
ratio.Divide(rRawMass)
ratio.Draw("P")
ratio.SetMarkerStyle(20)
ratio.GetYaxis().SetTitle("thr=0/thr=15")
ratio.GetYaxis().SetRangeUser(0.95, 1.05)
pad2.SetGridy(1)

canvas.Update()
canvas.SaveAs("thresholdComparison_rawMass.pdf")



 
