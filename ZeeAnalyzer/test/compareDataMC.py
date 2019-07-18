import ROOT 

ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

#no recovery, pierre tag, PU reweight
#mcFile=ROOT.TFile.Open("electronTree_PierreTag_MCPU_0GeV_noRecov_reweight.root", "READ")
#dataFile=ROOT.TFile.Open("electronTreeZEE_PierreTag_noChange.root", "READ")

#no recovery, pierre tag, data in 10_1_9 
#mcFile=ROOT.TFile.Open("", "READ")
#dataFile=ROOT.TFile.Open("electronTreeZEE_noChanges_pierre_1019.root", "READ")

#official MC, official Data (no pierre tag), CMSSW_10_2_X
#mcFile=ROOT.TFile.Open("electronTreeZEE_officialAODSIMFlat0to70PU.root", "READ")
#dataFile=ROOT.TFile.Open("electronTreeZEE_official_xtalInfo.root", "READ")

#official MC, data in 10_1_9
#mcFile=ROOT.TFile.Open("electronTreeZEE_officialAODSIMFlat0to70PU.root", "READ")
#dataFile=ROOT.TFile.Open("electronTreeZEE_noChanges_1019.root", "READ")

##recovery, pierre tag, sum8>0, xtalTrh>0.7, PU reweight
mcFile=ROOT.TFile.Open("electronTree_PierreTag_MCPU_0GeV_reweight.root", "READ")
dataFile=ROOT.TFile.Open("electronTreeZEE_PierreTag_sum8gt0_xtalInfo.root", "READ")

mcNtuple=mcFile.Get("ntupler/selected")
print mcNtuple.GetEntries()
dataNtuple=dataFile.Get("ntupler/selected")
print mcNtuple.GetEntries(), dataNtuple.GetEntries()

mc_e1R9=ROOT.TH1F("mc_e1R9", "mc_e1R9", 25, 0, 1)
data_e1R9=ROOT.TH1F("data_e1R9", "data_e1R9", 25, 0, 1)

mc_e2R9=ROOT.TH1F("mc_e2R9", "mc_e2R9", 25, 0, 1)
data_e2R9=ROOT.TH1F("data_e2R9", "data_e2R9", 25, 0, 1)

for evt in mcNtuple:
    ##if (evt.e1SeedIEta==-99 a
    if not (evt.e1SeedIEta==-99 or evt.e1SeedIPhi==-599) and evt.e1IsRecovered:
        mc_e1R9.Fill(evt.e1R9)
       
    if not (evt.e2SeedIEta==-99 or evt.e2SeedIPhi==-599) and evt.e2IsRecovered:
        mc_e2R9.Fill(evt.e2R9)

for evt in dataNtuple:
    if not (evt.e1SeedIEta==-99 or evt.e1SeedIPhi==-599) and evt.e1IsRecovered:
        data_e1R9.Fill(evt.e1R9)
    if not (evt.e2SeedIEta==-99 or evt.e2SeedIPhi==-599) and evt.e2IsRecovered:
        data_e2R9.Fill(evt.e2R9)


canvas=ROOT.TCanvas("c", "c", 600, 600)
print data_e1R9.Integral()
data_e1R9.Scale(1./data_e1R9.Integral())
data_e2R9.Scale(1./data_e2R9.Integral())
mc_e1R9.Scale(1./mc_e1R9.Integral())
mc_e2R9.Scale(1./mc_e2R9.Integral())

data_e1R9.Draw()
mc_e1R9.Draw("SAMES")
mc_e1R9.SetLineColor(2)
data_e1R9.GetXaxis().SetTitle(data_e1R9.GetTitle().replace('data_', ''))
mymax = max(mc_e1R9.GetBinContent(mc_e1R9.GetMaximumBin()), data_e1R9.GetBinContent(data_e1R9.GetMaximumBin()))
data_e1R9.GetYaxis().SetRangeUser(0,1.1*mymax)

legend=ROOT.TLegend(0.25, 0.8, 0.5, 0.6)
legend.AddEntry(data_e1R9, "data", "l")
legend.AddEntry(mc_e1R9, "mc", "l")

legend.Draw()

canvas.SaveAs("e1R9.png")
canvas.SaveAs("e1R9.root")

canvas.Clear()

data_e2R9.Draw()
mc_e2R9.Draw("SAMES")
mc_e2R9.SetLineColor(2)

mymax = max(mc_e2R9.GetBinContent(mc_e2R9.GetMaximumBin()), data_e2R9.GetBinContent(data_e2R9.GetMaximumBin()))
data_e2R9.GetYaxis().SetRangeUser(0,1.1*mymax)
data_e2R9.GetXaxis().SetTitle(data_e2R9.GetTitle().replace('data_', ''))
legend.Draw()


canvas.SaveAs("e2R9.png")
canvas.SaveAs("e2R9.root")
