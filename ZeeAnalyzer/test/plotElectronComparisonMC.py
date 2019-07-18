import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

infile=ROOT.TFile.Open("outfileEleComparisonDY.root", "READ")

histoList = [
'InvMass',
'InvMassRaw',
'e1R9',
'e2R9',
'e1RawEnergy',
'e2RawEnergy',
'e1Energy',
'e2Energy',
'e1SigmaIetaIeta',
'e2SigmaIetaIeta',
'e1RawGenRatio',
'e2RawGenRatio',
'e1EnGenRatio',
'e2EnGenRatio'
]

htitles=[
'Mass (GeV)',
'Raw Mass (GeV)',
'e1 R9', 
'e2 R9', 
'e1 raw Energy (GeV)',
'e2 raw Energy (GeV)', 
'e1 Energy (GeV)',
'e2 Energy (GeV)', 
'e1 #sigma i#eta i#eta', 
'e2 #sigma i#eta i#eta',
'e1 raw En / gen En',
'e2 raw En / gen En',
'e1 En / gen En',
'e2 En / gen En'
]

canvas=ROOT.TCanvas("c", "c", 600, 600)


for i,hname in enumerate(histoList):
    xmin=0.55
    if 'R9' in hname:
        xmin=0.45
    elif 'RawEnergy' in hname:
        xmin=0.65
    ymax=0.85
    xmax=xmin+0.2
    ymin=ymax-0.15


    legend = ROOT.TLegend(xmin, ymax, xmax, ymin)
    
    legend.Clear()

    mymax=0
    hP = infile.Get('p'+hname)
    hR = infile.Get('r'+hname)

    hP.Draw()
    hR.Draw('SAME')
    
    mymax=hP.GetBinContent(hP.GetMaximumBin())
    if mymax<hR.GetBinContent(hR.GetMaximumBin()): mymax=hR.GetBinContent(hR.GetMaximumBin())

    canvas.SetLogy(0)
    if 'SigmaIetaIeta' in hname: 
        canvas.SetLogy(1)
    hP.GetYaxis().SetRangeUser(0.1, 1.2*mymax)

    hP.GetXaxis().SetTitle(htitles[i])

    textxmin=0.15
    textymax=0.85
    if 'SigmaIetaIeta' in hname:
        textxmin=0.55
        textymax=0.65
    elif 'RawEnergy' in hname:
        textxmin=0.65
        textymax=0.65
    textxmax=textxmin+0.2
    textymin=textymax-0.1


    rounddigit = 2
    if 'SigmaIetaIeta' in hname:
        rounddigit=4

    rText=ROOT.TPaveText(textxmin, textymax, textxmax, textymin,"NDC")
    
    rText.AddText("After recovery:")
    rText.AddText("mean = "+str(round(hR.GetMean(),rounddigit)))
    rText.AddText("rms = "+str(round(hR.GetRMS(),rounddigit)));    
    rText.SetTextSize(rText.GetTextSize()*0.7);
    rText.SetTextAlign(11)
    rText.SetTextColor(2)    
    rText.SetFillColor(0)
    rText.Draw("SAME")

    pText=ROOT.TPaveText(textxmin, textymax-0.12, textxmax, textymin-0.12,"NDC")

    pText.AddText("As in production:")
    pText.AddText("mean = "+str(round(hP.GetMean(),rounddigit)))
    pText.AddText("rms = "+str(round(hP.GetRMS(),rounddigit)))    
    pText.SetTextSize(rText.GetTextSize()*0.7);
    pText.SetTextAlign(11)
    pText.SetTextColor(1)    
    pText.SetFillColor(0)
    pText.Draw("SAME")


    hP.SetLineColor(1)
    hP.SetLineWidth(2)
    hR.SetLineWidth(2)
    hR.SetLineColor(2)

    
    legend.AddEntry(hP, 'production', 'l')
    legend.AddEntry(hR, 'recovery', 'l')

    legend.Draw()
    

    canvas.Update()
    canvas.SaveAs('plotsDY/'+hname+'.png')
    canvas.SaveAs('plotsDY/'+hname+'.pdf')


infile.Close()
