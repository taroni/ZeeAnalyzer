import ROOT
import time

eventlist = [evt.rstrip('\n') for evt in open("events.txt")]
recovFile = ROOT.TFile.Open("electronTree_15GeV.root", "READ")
prodFile = ROOT.TFile.Open("electronTree_PierreTag.root", "READ")

rTree = recovFile.Get("ntupler/selected")
pTree = prodFile.Get("ntupler/selected")


outfile= ROOT.TFile.Open("outfileEleComparison_noThr_15GeV.root", "RECREATE")
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

count=0
for row in rTree:
    count+=1
    #print 'processing electron # :', count
    if count %1000 == 1:
        now=time.strftime("%H:%M:%S", time.localtime())
        print 'processing the %sst electron at %s' %(str(count), now)
    if row.e1Charge*row.e2Charge!=-1: continue
    if (row.e1SeedIPhi==-599.) : continue
    if (row.e2SeedIPhi==-599.) : continue
    if (row.e1SeedIEta==-99.) : continue
    if (row.e2SeedIEta==-99.) : continue

    myevt='%s:%s:%s' %(str(row.runNumber), str(row.lumiBlock), str(row.eventNumber))
    if myevt in eventlist:
                print 'event %s:%s:%s' %(str(row.runNumber), str(row.lumiBlock), str(row.eventNumber))
                print 'mass (raw SC)',  row.invMass_rawSC 
                print 'raw energy e1', row.e1RawEnergy,'e2', row.e2RawEnergy
                print 'energy e1', row.e1Energy, 'e2', row.e2Energy
                print 'charge e1:', row.e1Charge, 'e2', row.e2Charge
                print 'phiSC e1:', row.e1PhiSC, 'e2', row.e2PhiSC
                print 'etaSC e1:',row.e1EtaSC,  'e2', row.e2EtaSC
                print 'phi e1:', row.e1Phi,  'e2', row.e2Phi
                print 'eta e1:',row.e1Eta, 'e2', row.e2Eta
    for entry in pTree:
        if entry.e1Charge*entry.e2Charge!=-1 : continue
        if not row.runNumber==entry.runNumber: continue
        if not row.lumiBlock==entry.lumiBlock: continue
        if not row.eventNumber==entry.eventNumber: continue
        if (entry.e1SeedIPhi==-599.) : continue
        if (entry.e2SeedIPhi==-599.) : continue
        if (entry.e1SeedIEta==-99.) : continue
        if (entry.e2SeedIEta==-99.) : continue

             
        if myevt in eventlist:
            print 'event %s:%s:%s' %(str(row.runNumber), str(row.lumiBlock), str(row.eventNumber))
            print 'mass (raw SC)',  row.invMass_rawSC , entry.invMass_rawSC 
            print 'raw energy e1', row.e1RawEnergy, entry.e1RawEnergy, 'e2', row.e2RawEnergy, entry.e2RawEnergy
            print 'energy e1', row.e1Energy, entry.e1Energy, 'e2', row.e2Energy, entry.e2Energy
            print 'charge e1:', row.e1Charge,  entry.e1Charge, 'e2', row.e2Charge, entry.e2Charge
            print 'phiSC e1:', row.e1PhiSC, entry.e1PhiSC, 'e2', row.e2PhiSC, entry.e2PhiSC
            print 'etaSC e1:',row.e1EtaSC, entry.e1EtaSC, 'e2', row.e2EtaSC, entry.e2EtaSC
            print 'iphiSeed e1:', row.e1SeedIPhi, entry.e1SeedIPhi, 'e2', row.e2SeedIPhi, entry.e2SeedIPhi
            print 'ietaSeed e1:' ,row.e1SeedIEta, entry.e1SeedIEta, 'e2', row.e2SeedIEta, entry.e2SeedIEta
            print 'phi e1:', row.e1Phi, entry.e1Phi, 'e2', row.e2Phi, entry.e2Phi
            print 'eta e1:',row.e1Eta, entry.e1Eta, 'e2', row.e2Eta, entry.e2Eta
            
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

print 'comparison ended'

outfile.Close()
