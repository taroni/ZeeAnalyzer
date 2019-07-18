import ROOT
import ast

isMC=False
#recovFile = ROOT.TFile.Open("electronTree_PierreTag_MCPU_0GeV_reweight.root", "READ")
#prodFile = ROOT.TFile.Open("electronTree_PierreTag_MCPU_0GeV_noRecov_reweight.root", "READ")
#prodFile = ROOT.TFile.Open("electronTreeZEE_PierreTag_noChange.root", "READ")
prodFile = ROOT.TFile.Open("electronTreeZEE_NoChange_sum8gt20_xtalInfo.root", "READ")
recovFile = ROOT.TFile.Open("electronTreeZEE_PierreTag_sum8gt20_xtalInfo.root", "READ")

rTree = recovFile.Get("ntupler/selected")
pTree = prodFile.Get("ntupler/selected")


outfile= ROOT.TFile.Open("outfileEleInCluster_data.root", "RECREATE")

rInvMassRaw=ROOT.TH1F("rInvMassRaw", "rInvMassRaw", 50, 0, 180)
pInvMassRaw=ROOT.TH1F("pInvMassRaw", "pInvMassRaw", 50, 0, 180)

h0=ROOT.TH1F("singleClustXtalEn_recovProdDiff","singleClustXtalEn_recovProdDiff", 200, -100, 100)
h1=ROOT.TH1F("twoClustXtalEn_recovProdDiff","twoClustXtalEn_recovProdDiff", 200, -100, 100)

h0Fr=ROOT.TH1F("singleClustXtalEnFr_recovProdDiff","singleClustXtalEnFr_recovProdDiff", 200, -5, 5)
h1Fr=ROOT.TH1F("twoClustXtalEnFr_recovProdDiff","twoClustXtalEnFr_recovProdDiff", 200, -5, 5)

hSingleClustEnRecov=ROOT.TH1F("singleClustXtalEnRecov", "singleClustXtalEnRecov", 200, 0, 200)
hTwoClustEnRecov=ROOT.TH1F("twoClustXtalEnRecov", "twoClustXtalEnRecov",200, 0, 200)
hSingleClustEnProd=ROOT.TH1F("singleClustXtalEnProd", "singleClustXtalEnProd", 200, 0, 200)
hTwoClustEnProd=ROOT.TH1F("twoClustXtalEnProd", "twoClustXtalEnProd",200, 0, 200)

singleClustDiffvsEnProd=ROOT.TH2F("singleClustDiffvsEnProd", "singleClustDiffvsEnProd", 200, 0, 200, 200, -100, 100)
singleClustDiffvsEnRecov=ROOT.TH2F("singleClustDiffvsEnRecov", "singleClustDiffvsEnRecov", 200, 0, 200, 200, -100, 100)
twoClustDiffvsEnProd=ROOT.TH2F("twoClustDiffvsEnProd", "twoClustDiffvsEnProd", 200, 0, 200, 200, -100, 100)
twoClustDiffvsEnRecov=ROOT.TH2F("twoClustDiffvsEnRecov", "twoClustDiffvsEnRecov", 200, 0, 200, 200, -100, 100)
singleClustDiffFrvsEnProd=ROOT.TH2F("singleClustDiffFrvsEnProd", "singleClustDiffFrvsEnProd", 200, 0, 200, 200, -5, 5)
singleClustDiffFrvsEnRecov=ROOT.TH2F("singleClustDiffFrvsEnRecov", "singleClustDiffFrvsEnRecov", 200, 0, 200, 200, -5, 5)
twoClustDiffFrvsEnProd=ROOT.TH2F("twoClustDiffFrvsEnProd", "twoClustDiffFrvsEnProd", 200, 0, 200, 200, -5, 5)
twoClustDiffFrvsEnRecov=ROOT.TH2F("twoClustDiffFrvsEnRecov", "twoClustDiffFrvsEnRecov", 200, 0, 200, 200, -5, 5)


s = open('Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt', 'r').read()
mydict=ast.literal_eval(s)

count=0
print 'total number of entries:', rTree.GetEntries() 
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
            

            if len(row.xtalEn)==1:
                #if len(entry.xtalEn)>1:
                #    print 'recovered', 0, row.xtalRawId[0], row.xtalEn[0]
                #    for jxtal in xrange(len(entry.xtalEn)):
                #         print 'official', jxtal, entry.xtalRawId[jxtal], entry.xtalEn[jxtal]
                if row.xtalRawId[0] == entry.xtalRawId[0]:
                    h0.Fill(row.xtalEn[0]-entry.xtalEn[0])
                    h0Fr.Fill((row.xtalEn[0]-entry.xtalEn[0])/entry.xtalEn[0])
                    hSingleClustEnRecov.Fill(row.xtalEn[0])
                    hSingleClustEnProd.Fill(entry.xtalEn[0])
                    
                    singleClustDiffvsEnProd.Fill(entry.xtalEn[0], row.xtalEn[0]-entry.xtalEn[0])
                    singleClustDiffvsEnRecov.Fill(row.xtalEn[0],  row.xtalEn[0]-entry.xtalEn[0])
                    singleClustDiffFrvsEnProd.Fill(entry.xtalEn[0], (row.xtalEn[0]-entry.xtalEn[0])/entry.xtalEn[0])
                    singleClustDiffFrvsEnRecov.Fill(row.xtalEn[0],  (row.xtalEn[0]-entry.xtalEn[0])/row.xtalEn[0])

            if len(row.xtalEn)>1:
                #print len(row.xtalEn), len(entry.xtalEn)
                #if len(row.xtalEn)!=len(entry.xtalEn):
                #     for ixtal in xrange(len(row.xtalEn)):
                #         print 'recovered', ixtal, row.xtalRawId[ixtal], row.xtalEn[ixtal]
                #     for jxtal in xrange(len(entry.xtalEn)):
                #         print 'official', jxtal, entry.xtalRawId[jxtal], entry.xtalEn[jxtal]
                for ixtal in xrange(len(row.xtalEn)):
                    for jxtal in xrange(len(entry.xtalEn)):
                        if row.xtalRawId[ixtal] == entry.xtalRawId[jxtal]:
                            h1.Fill(row.xtalEn[ixtal]-entry.xtalEn[jxtal])
                            h1Fr.Fill((row.xtalEn[ixtal]-entry.xtalEn[jxtal])/entry.xtalEn[jxtal])                            
                            hTwoClustEnRecov.Fill(row.xtalEn[ixtal])
                            hTwoClustEnProd.Fill(entry.xtalEn[jxtal])
                    
                            twoClustDiffvsEnProd.Fill(entry.xtalEn[jxtal], row.xtalEn[ixtal]-entry.xtalEn[jxtal])
                            twoClustDiffvsEnRecov.Fill(row.xtalEn[ixtal],  row.xtalEn[ixtal]-entry.xtalEn[jxtal])
                            twoClustDiffFrvsEnProd.Fill(entry.xtalEn[0], (row.xtalEn[0]-entry.xtalEn[0])/entry.xtalEn[0])
                            twoClustDiffFrvsEnRecov.Fill(row.xtalEn[0],  (row.xtalEn[0]-entry.xtalEn[0])/row.xtalEn[0])


print 'all event processed'      
outfile.cd()
h0.Write()
h1.Write()
h0Fr.Write()
h1Fr.Write()
hSingleClustEnRecov.Write()
hTwoClustEnRecov.Write()
hSingleClustEnProd.Write()
hTwoClustEnProd.Write()
singleClustDiffvsEnProd.Write()
singleClustDiffvsEnRecov.Write()
twoClustDiffvsEnProd.Write()
twoClustDiffvsEnRecov.Write()
singleClustDiffFrvsEnProd.Write()
singleClustDiffFrvsEnRecov.Write()
twoClustDiffFrvsEnProd.Write()
twoClustDiffFrvsEnRecov.Write()

outfile.Close()
print 'comparison completed'
