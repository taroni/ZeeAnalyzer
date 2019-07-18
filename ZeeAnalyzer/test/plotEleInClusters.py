import ROOT

ROOT.gROOT.SetBatch(1)

file0=ROOT.TFile.Open("outfileEleInCluster_data.root","READ")
#file0=ROOT.TFile.Open("output_Neighbours_noChange.root","READ")

canvas= ROOT.TCanvas("c","c", 600, 600)

canvas.Draw()

ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetHistLineWidth(2)
ROOT.gROOT.ForceStyle()

singleClustXtalEn_recovProdDiff=file0.Get("singleClustXtalEn_recovProdDiff")
twoClustXtalEn_recovProdDiff=file0.Get("twoClustXtalEn_recovProdDiff")

singleClustXtalEnFr_recovProdDiff=file0.Get("singleClustXtalEnFr_recovProdDiff")
twoClustXtalEnFr_recovProdDiff=file0.Get("twoClustXtalEnFr_recovProdDiff")

singleClustDiffvsEnProd=file0.Get("singleClustDiffvsEnProd")
singleClustDiffvsEnRecov=file0.Get("singleClustDiffvsEnRecov")

twoClustDiffvsEnProd=file0.Get("twoClustDiffvsEnProd")
twoClustDiffvsEnRecov=file0.Get("twoClustDiffvsEnRecov")

singleClustXtalEn_recovProdDiff.Rebin(5)
twoClustXtalEn_recovProdDiff.Rebin(5)
singleClustDiffvsEnProd.Rebin2D(2,2)
singleClustDiffvsEnRecov.Rebin2D(2,2)
twoClustDiffvsEnProd.Rebin2D(2,2)
twoClustDiffvsEnRecov.Rebin2D(2,2)

leg=ROOT.TLegend(0.55, 0.85, 0.85, 0.7)

canvas.SetRightMargin(0.15)
canvas.SetLeftMargin(0.13)
singleClustXtalEn_recovProdDiff.Scale(1./singleClustXtalEn_recovProdDiff.Integral())
singleClustXtalEn_recovProdDiff.Draw("HIST")
twoClustXtalEn_recovProdDiff.Scale(1./twoClustXtalEn_recovProdDiff.Integral())
twoClustXtalEn_recovProdDiff.Draw("HISTSAMES")
twoClustXtalEn_recovProdDiff.SetLineColor(2)

mymax=1.2*max(singleClustXtalEn_recovProdDiff.GetBinContent(singleClustXtalEn_recovProdDiff.GetMaximumBin()), twoClustXtalEn_recovProdDiff.GetBinContent(twoClustXtalEn_recovProdDiff.GetMaximumBin()))
singleClustXtalEn_recovProdDiff.GetXaxis().SetTitle("recovered - production crystal energy [GeV]")
singleClustXtalEn_recovProdDiff.GetYaxis().SetRangeUser(0, mymax)
leg.AddEntry(singleClustXtalEn_recovProdDiff, "one cluster", "l")
leg.AddEntry(twoClustXtalEn_recovProdDiff, "two (or more) clusters", "l")
leg.Draw()
canvas.Update()
canvas.SaveAs("plots/xtalEn_recovProdDiff.png")
canvas.SaveAs("plots/xtalEn_recovProdDiff.pdf")
canvas.Clear()
singleClustXtalEnFr_recovProdDiff.Scale(1./singleClustXtalEnFr_recovProdDiff.Integral())
singleClustXtalEnFr_recovProdDiff.Draw("HIST")
twoClustXtalEnFr_recovProdDiff.Scale(1./twoClustXtalEnFr_recovProdDiff.Integral())
twoClustXtalEnFr_recovProdDiff.Draw("HISTSAMES")
twoClustXtalEnFr_recovProdDiff.SetLineColor(2)
mymax=1.2*max(singleClustXtalEnFr_recovProdDiff.GetBinContent(singleClustXtalEnFr_recovProdDiff.GetMaximumBin()), twoClustXtalEnFr_recovProdDiff.GetBinContent(twoClustXtalEnFr_recovProdDiff.GetMaximumBin()))
singleClustXtalEnFr_recovProdDiff.GetXaxis().SetTitle("recovered - production crystal energy [GeV]")
singleClustXtalEnFr_recovProdDiff.GetYaxis().SetRangeUser(0, mymax)
leg.Draw()
canvas.Update()
canvas.SaveAs("plots/xtalEnFr_recovProdDiff.png")
canvas.SaveAs("plots/xtalEnFr_recovProdDiff.pdf")

canvas.Clear()
singleClustDiffvsEnProd.Draw("COLZ")
singleClustDiffvsEnProd.GetYaxis().SetTitle("recovered - production crystal energy [GeV]")
singleClustDiffvsEnProd.GetXaxis().SetTitle("production xtal energy [GeV]")
canvas.SaveAs("plots/singleClust_xtalEn_recovProdDiffvsProdEn.png")
canvas.SaveAs("plots/singleClust_xtalEn_recovProdDiffvsProdEn.pdf")

canvas.Clear()
twoClustDiffvsEnProd.Draw("COLZ")
twoClustDiffvsEnProd.GetYaxis().SetTitle("recovered - production crystal energy [GeV]")
twoClustDiffvsEnProd.GetXaxis().SetTitle("production xtal energy [GeV]")
canvas.SaveAs("plots/twoClust_xtalEn_recovProdDiffvsProdEn.png")
canvas.SaveAs("plots/twoClust_xtalEn_recovProdDiffvsProdEn.pdf")


file0.Close()


