import ROOT
import ast

isMC=False
#recovFile = ROOT.TFile.Open("electronTree_PierreTag_MCPU_0GeV_reweight.root", "READ")
#prodFile = ROOT.TFile.Open("electronTree_PierreTag_MCPU_0GeV_noRecov_reweight.root", "READ")
#prodFile = ROOT.TFile.Open("electronTreeZEE_PierreTag_noChange.root", "READ")
prodFile = ROOT.TFile.Open("electronTree_PierreTag_0_190412.root", "READ")
recovFile = ROOT.TFile.Open("electronTreeZEE_PierreTag_sum8gt0_xtalInfo.root", "READ")

rTree = recovFile.Get("ntupler/selected")
pTree = prodFile.Get("ntupler/selected")


outfile= ROOT.TFile.Open("outfileEleComparison_data.root", "RECREATE")
##rTree.SetBranchAddress("runNumber", rRun)
##rTree.SetBranchAddress("lumiBlock", rLumi)
##rTree.SetBranchAddress("eventNumber", rEvt)
##pTree.SetBranchAddress("runNumber", pRun)
##pTree.SetBranchAddress("lumiBlock", pLumi)
##pTree.SetBranchAddress("eventNumber", pEvt)

rInvMass=ROOT.TH1F("rInvMass", "rInvMass", 50, 0, 180)
pInvMass=ROOT.TH1F("pInvMass", "pInvMass", 50, 0, 180)
rInvMassRaw=ROOT.TH1F("rInvMassRaw", "rInvMassRaw", 50, 0, 180)
pInvMassRaw=ROOT.TH1F("pInvMassRaw", "pInvMassRaw", 50, 0, 180)

re1R9=ROOT.TH1F("re1R9", "re1R9", 24, 0, 1.2)
pe1R9=ROOT.TH1F("pe1R9", "pe1R9", 24, 0, 1.2)
re2R9=ROOT.TH1F("re2R9", "re2R9", 24, 0, 1.2)
pe2R9=ROOT.TH1F("pe2R9", "pe2R9", 24, 0, 1.2)

re1RawEnergy=ROOT.TH1F("re1RawEnergy", "re1RawEnergy", 40, 0, 200)
pe1RawEnergy=ROOT.TH1F("pe1RawEnergy", "pe1RawEnergy", 40, 0, 200)
re2RawEnergy=ROOT.TH1F("re2RawEnergy", "re2RawEnergy", 40, 0, 200)
pe2RawEnergy=ROOT.TH1F("pe2RawEnergy", "pe2RawEnergy", 40, 0, 200)
re1Energy=ROOT.TH1F("re1Energy", "re1Energy", 40, 0, 200)
pe1Energy=ROOT.TH1F("pe1Energy", "pe1Energy", 40, 0, 200)
re2Energy=ROOT.TH1F("re2Energy", "re2Energy", 40, 0, 200)
pe2Energy=ROOT.TH1F("pe2Energy", "pe2Energy", 40, 0, 200)

re1SigmaIetaIeta=ROOT.TH1F("re1SigmaIetaIeta","re1SigmaIetaIeta", 100, 0, 0.1)
re2SigmaIetaIeta=ROOT.TH1F("re2SigmaIetaIeta","re2SigmaIetaIeta", 100, 0, 0.1)
pe1SigmaIetaIeta=ROOT.TH1F("pe1SigmaIetaIeta","pe1SigmaIetaIeta", 100, 0, 0.1)
pe2SigmaIetaIeta=ROOT.TH1F("pe2SigmaIetaIeta","pe2SigmaIetaIeta", 100, 0, 0.1)

s = open('Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt', 'r').read()
mydict=ast.literal_eval(s)

count=0
for row in rTree:

    if count%100==0: print 'the %ith entry has been analysed' %count
    count+=1
    if not isMC:
        if not  str(row.runNumber) in mydict: continue
        lumilist=mydict[str(row.runNumber)]
        lumiCheck=False
        for lumiInt in lumilist:
            start, end = lumiInt
            if row.lumiBlock in range (start, end+1):
                lumiCheck=True
                break
        if lumiCheck==False: continue

    if row.e1Charge*row.e2Charge!=-1: continue
    if (row.e1SeedIPhi==-599. and row.e2SeedIPhi==-599.) : continue
    if (row.e1SeedIEta==-99. and row.e2SeedIEta==-99.) : continue
    ##if (row.e1SeedIPhi==-599.) : continue
    ##if (row.e2SeedIPhi==-599.) : continue
    ##if (row.e1SeedIEta==-99.) : continue
    ##if (row.e2SeedIEta==-99.) : continue

    for entry in pTree:

        if not row.runNumber==entry.runNumber: continue
        if not row.lumiBlock==entry.lumiBlock: continue
        if not row.eventNumber==entry.eventNumber: continue
        if entry.e1Charge*entry.e2Charge!=-1 : continue
        if (entry.e1SeedIPhi==-599. and entry.e2SeedIPhi==-599.) : continue
        if (entry.e1SeedIEta==-99. and entry.e2SeedIEta==-99.) : continue

        if (entry.e1SeedIPhi==-599.) : continue
        if (entry.e2SeedIPhi==-599.) : continue
        if (entry.e1SeedIEta==-99.) : continue
        if (entry.e2SeedIEta==-99.) : continue

        if  ( abs(row.e1SeedIPhi-entry.e1SeedIPhi)>1 and abs(row.e2SeedIPhi-entry.e2SeedIPhi)>1 and  abs(row.e1SeedIEta-entry.e1SeedIEta)>1 and abs(row.e2SeedIEta-entry.e2SeedIEta)>1) or (  abs(row.e1SeedIPhi-entry.e2SeedIPhi)>1 and abs(row.e2SeedIPhi-entry.e1SeedIPhi)>1 and  abs(row.e1SeedIEta-entry.e2SeedIEta)>1 and abs(row.e2SeedIEta-entry.e1SeedIEta)>1) :
            if abs(row.e1EtaSC-entry.e1EtaSC)>0.1 : continue
            if abs(row.e2EtaSC-entry.e2EtaSC)>0.1 : continue
            if abs(row.e1PhiSC-entry.e1PhiSC)>0.1 : continue
            if abs(row.e2PhiSC-entry.e2PhiSC)>0.1 : continue

            if not (row.e1IsRecovered or row.e2IsRecovered): continue 
            

            rInvMass.Fill(row.invMass)
            rInvMassRaw.Fill(row.invMass_rawSC)
            
            pInvMass.Fill(entry.invMass)
            pInvMassRaw.Fill(entry.invMass_rawSC)
            
            if (row.e1IsRecovered):
                re1R9.Fill(row.e1R9)
                re1RawEnergy.Fill(row.e1RawEnergy)
                re1Energy.Fill(row.e1Energy)
                re1SigmaIetaIeta.Fill(row.e1SigmaIetaIeta)
                pe1R9.Fill(entry.e1R9)
                pe1RawEnergy.Fill(entry.e1RawEnergy)
                pe1Energy.Fill(entry.e1Energy)
                pe1SigmaIetaIeta.Fill(entry.e1SigmaIetaIeta)
            
            elif(row.e2IsRecovered):
                re2R9.Fill(row.e2R9)
                re2RawEnergy.Fill(row.e2RawEnergy)
                re2Energy.Fill(row.e2Energy)
                re2SigmaIetaIeta.Fill(row.e2SigmaIetaIeta)
                pe2R9.Fill(entry.e2R9)
                pe2RawEnergy.Fill(entry.e2RawEnergy)
                pe2Energy.Fill(entry.e2Energy)
                pe2SigmaIetaIeta.Fill(entry.e2SigmaIetaIeta)
                
print 'all event processed'      
outfile.cd()
rInvMass.Write()
pInvMass.Write()
rInvMassRaw.Write()
pInvMassRaw.Write()

re1R9.Write()
pe1R9.Write()
re2R9.Write()
pe2R9.Write()

re1RawEnergy.Write()
pe1RawEnergy.Write()
re2RawEnergy.Write()
pe2RawEnergy.Write()
re1Energy.Write()
pe1Energy.Write()
re2Energy.Write()
pe2Energy.Write()

re1SigmaIetaIeta.Write()
re2SigmaIetaIeta.Write()
pe1SigmaIetaIeta.Write()
pe2SigmaIetaIeta.Write()

outfile.Close()
print 'comparison completed'
