import ROOT

recovFile = ROOT.TFile.Open("electronTreeDY.root", "READ")
prodFile = ROOT.TFile.Open("electronTreeDYnoChanges.root", "READ")

rTree = recovFile.Get("ntupler/selected")
pTree = prodFile.Get("ntupler/selected")


outfile= ROOT.TFile.Open("outfileEleComparisonDY.root", "RECREATE")
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

re1GenEnergy=ROOT.TH1F("re1GenEnergy", "re1GenEnergy", 40, 0, 200)
pe1GenEnergy=ROOT.TH1F("pe1GenEnergy", "pe1GenEnergy", 40, 0, 200)
re2GenEnergy=ROOT.TH1F("re2GenEnergy", "re2GenEnergy", 40, 0, 200)
pe2GenEnergy=ROOT.TH1F("pe2GenEnergy", "pe2GenEnergy", 40, 0, 200)

re1RawGenRatio=ROOT.TH1F("re1RawGenRatio", "re1RawGenRatio", 40, 0, 2)
pe1RawGenRatio=ROOT.TH1F("pe1RawGenRatio", "pe1RawGenRatio", 40, 0, 2)
re2RawGenRatio=ROOT.TH1F("re2RawGenRatio", "re2RawGenRatio", 40, 0, 2)
pe2RawGenRatio=ROOT.TH1F("pe2RawGenRatio", "pe2RawGenRatio", 40, 0, 2)
re1EnGenRatio=ROOT.TH1F("re1EnGenRatio", "re1EnGenRatio", 40, 0, 2)
pe1EnGenRatio=ROOT.TH1F("pe1EnGenRatio", "pe1EnGenRatio", 40, 0, 2)
re2EnGenRatio=ROOT.TH1F("re2EnGenRatio", "re2EnGenRatio", 40, 0, 2)
pe2EnGenRatio=ROOT.TH1F("pe2EnGenRatio", "pe2EnGenRatio", 40, 0, 2)


re1SigmaIetaIeta=ROOT.TH1F("re1SigmaIetaIeta","re1SigmaIetaIeta", 100, 0, 0.1)
re2SigmaIetaIeta=ROOT.TH1F("re2SigmaIetaIeta","re2SigmaIetaIeta", 100, 0, 0.1)
pe1SigmaIetaIeta=ROOT.TH1F("pe1SigmaIetaIeta","pe1SigmaIetaIeta", 100, 0, 0.1)
pe2SigmaIetaIeta=ROOT.TH1F("pe2SigmaIetaIeta","pe2SigmaIetaIeta", 100, 0, 0.1)

count=0
for row in rTree:
    count+=1
    if count % 1000 == 1 : print 'processing %sst electron' % str(count)
    if row.e1Charge*row.e2Charge!=-1: continue
    if (row.e1SeedIPhi==-599.) : continue
    if (row.e2SeedIPhi==-599.) : continue
    if (row.e1SeedIEta==-99.) : continue
    if (row.e2SeedIEta==-99.) : continue

    for entry in pTree:
        if entry.e1Charge*entry.e2Charge!=-1 : continue
        if not row.runNumber==entry.runNumber: continue
        if not row.lumiBlock==entry.lumiBlock: continue
        if not row.eventNumber==entry.eventNumber: continue
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
            
            #print row.runNumber, row.lumiBlock, row.eventNumber, entry.runNumber, entry.lumiBlock, entry.eventNumber, row.invMass, entry.invMass, row.e1Charge, row.e2Charge, entry.e1Charge, entry.e2Charge, row.e1PhiSC, entry.e1PhiSC, row.e2PhiSC, entry.e2PhiSC

            rInvMass.Fill(row.invMass)
            rInvMassRaw.Fill(row.invMass_rawSC)
            
            pInvMass.Fill(entry.invMass)
            pInvMassRaw.Fill(entry.invMass_rawSC)
 
            if (row.e1IsRecovered):
                re1R9.Fill(row.e1R9)
                re1RawEnergy.Fill(row.e1RawEnergy)
                re1Energy.Fill(row.e1Energy)
                re1GenEnergy.Fill(row.e1GenEnergy)
                re1RawGenRatio.Fill(row.e1RawEnergy/row.e1GenEnergy)
                re1EnGenRatio.Fill(row.e1Energy/row.e1GenEnergy)
                re1SigmaIetaIeta.Fill(row.e1SigmaIetaIeta)
                pe1R9.Fill(entry.e1R9)
                pe1RawEnergy.Fill(entry.e1RawEnergy)
                pe1Energy.Fill(entry.e1Energy)
                pe1GenEnergy.Fill(entry.e1GenEnergy)
                pe1RawGenRatio.Fill(entry.e1RawEnergy/entry.e1GenEnergy)
                pe1EnGenRatio.Fill(entry.e1Energy/entry.e1GenEnergy)
                pe1SigmaIetaIeta.Fill(entry.e1SigmaIetaIeta)

            elif(row.e2IsRecovered):
                re2R9.Fill(row.e2R9)
                re2RawEnergy.Fill(row.e2RawEnergy)
                re2Energy.Fill(row.e2Energy)
                re2GenEnergy.Fill(row.e2GenEnergy)
                re2RawGenRatio.Fill(row.e2RawEnergy/row.e2GenEnergy)
                re2EnGenRatio.Fill(row.e2Energy/row.e2GenEnergy)
                re2SigmaIetaIeta.Fill(row.e2SigmaIetaIeta)
                pe2R9.Fill(entry.e2R9)
                pe2RawEnergy.Fill(entry.e2RawEnergy)
                pe2Energy.Fill(entry.e2Energy)
                pe2GenEnergy.Fill(entry.e2GenEnergy)
                pe2RawGenRatio.Fill(entry.e2RawEnergy/entry.e2GenEnergy)
                pe2EnGenRatio.Fill(entry.e2Energy/entry.e2GenEnergy)
                pe2SigmaIetaIeta.Fill(entry.e2SigmaIetaIeta)

            break
                
            
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

pe2GenEnergy.Write()
pe2RawGenRatio.Write()
pe2EnGenRatio.Write()
re2GenEnergy.Write()
re2RawGenRatio.Write()
re2EnGenRatio.Write()
pe1GenEnergy.Write()
pe1RawGenRatio.Write()
pe1EnGenRatio.Write()
re1GenEnergy.Write()
re1RawGenRatio.Write()
re1EnGenRatio.Write()


print 'comparison ended'
outfile.Close()
