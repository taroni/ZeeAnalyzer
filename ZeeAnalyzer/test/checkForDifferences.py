import ROOT
import uproot
import numpy as np
import pandas as pd
import array as arr

import argparse

parser = argparse.ArgumentParser(description='createHistograms')
parser.add_argument('--createHistograms', dest='histos', action='store', default=False, help='boolean if saving histos', required=False)

args = parser.parse_args()
createHistos=bool(args.histos)

tree0 =uproot.open("electronTree_PierreTag_MCPU_0GeV_reweight.root")['ntupler']['selected']
tree_noCh=uproot.open("electronTree_PierreTag_MCPU_0GeV_noRecov_reweight.root")['ntupler']['selected']
#tree0 =uproot.open("electronTreeZEE_PierreTag_sum8gt0_xtalInfo.root")['ntupler']['selected']
#tree_noCh=uproot.open("electronTreeZEE_PierreTag_noChange.root")['ntupler']['selected']

#outfile= ROOT.TFile.Open("outfileEleComparisonDY_0_15GeV.root", "RECREATE")

indata0=tree0.pandas.df(['runNumber', 'lumiBlock', 'eventNumber','e1Charge','e2Charge','e1Eta', 'e1Phi', 'e1R9', 'e2Eta', 'e2Phi', 'e2R9', 'e1SeedIEta', 'e1SeedIPhi', 'e2SeedIEta', 'e2SeedIPhi', 'e1SigmaIetaIeta', 'e2SigmaIetaIeta', 'e1EtaSC', 'e2EtaSC', 'e1PhiSC', 'e2PhiSC',   'e1Energy', 'e1RawEnergy', 'e2Energy', 'e2RawEnergy',  'invMass', 'invMass_rawSC', 'e1GenEnergy', 'e2GenEnergy', 'e1IsDead', 'e1IsRecovered', 'e2IsDead', 'e2IsRecovered'
])
print list(indata0.columns.values), indata0.shape[0]

indata_noCh=tree_noCh.pandas.df(['runNumber', 'lumiBlock', 'eventNumber', 'e1Charge', 'e2Charge', 'e1Eta', 'e1Phi', 'e1R9', 'e2Eta', 'e2Phi', 'e2R9', 'e1SeedIEta', 'e1SeedIPhi', 'e2SeedIEta', 'e2SeedIPhi', 'e1SigmaIetaIeta', 'e2SigmaIetaIeta', 'e1EtaSC', 'e2EtaSC', 'e1PhiSC', 'e2PhiSC', 'e1Energy', 'e1RawEnergy', 'e2Energy', 'e2RawEnergy', 'invMass', 'invMass_rawSC', 'e1GenEnergy', 'e2GenEnergy', 'e1IsDead', 'e1IsRecovered', 'e2IsDead', 'e2IsRecovered'
])


indata_noCh=indata_noCh[indata_noCh.e1Charge*indata_noCh.e2Charge==-1]
indata0=indata0[indata0.e1Charge*indata0.e2Charge==-1]

#indata_noCh=indata_noCh[(indata_noCh.invMass_rawSC>60)  &(indata_noCh.invMass_rawSC<120)]
#indata0=indata0[(indata0.invMass_rawSC>60)  &(indata0.invMass_rawSC<120)]

#indata_noCh=indata_noCh[indata_noCh.e1SeedIEta!=-99]
#indata_noCh=indata_noCh[indata_noCh.e1SeedIPhi!=-599]
#indata0=indata0[indata0.e1SeedIEta!=-99]
#indata0=indata0[indata0.e1SeedIPhi!=-599]

#indata_noCh=indata_noCh[(indata_noCh.e1IsRecovered==1) | (indata_noCh.e2IsRecovered==1)]
#indata0=indata0[(indata0.e1IsRecovered==1) | (indata0.e2IsRecovered==1)]
print list(indata0.columns.values), indata0.shape[0]
print list(indata_noCh.columns.values), indata_noCh.shape[0]

evtlist0 = set(["%i:%i:%i" %(indata0.iloc[x,:].runNumber, indata0.iloc[x,:].lumiBlock, indata0.iloc[x,:].eventNumber) for x in range(0, len(indata0))])
evtlist = set(["%i:%i:%i" %(indata_noCh.iloc[x,:].runNumber, indata_noCh.iloc[x,:].lumiBlock, indata_noCh.iloc[x,:].eventNumber) for x in range(0, len(indata_noCh))])
print len(evtlist0), len(evtlist)

#df=pd.merge(indata_noCh, indata0, how='inner',on=["runNumber", "lumiBlock", "eventNumber", 'e1Charge', 'e2Charge',  'e1SeedIEta', 'e1SeedIPhi', 'e2SeedIEta', 'e2SeedIPhi','e1GenEnergy', 'e2GenEnergy', 'e1IsRecovered', 'e2IsRecovered'])#, indicator=True) 
#print indata_noCh.head()
print indata0.head(20)
#removing gen energies when using real data
df=pd.merge(indata_noCh, indata0, how='inner',on=["runNumber", "lumiBlock", "eventNumber"], indicator=True) 


#df=df[(abs(df.e1PhiSC_x-df.e1PhiSC_y)<0.5) & (abs(df.e2PhiSC_x-df.e2PhiSC_y)<0.5) & (abs(df.e1EtaSC_x-df.e1EtaSC_y)<0.5) & (abs(df.e2EtaSC_x-df.e2EtaSC_y)<0.5)]
#print df

diffE1 = df.e1Energy_x-df.e1Energy_y
diffE2 = df.e2Energy_x-df.e2Energy_y
diffRawE1 = df.e1RawEnergy_x-df.e1RawEnergy_y
diffRawE2 = df.e2RawEnergy_x-df.e2RawEnergy_y
diffMass = df.invMass_x-df.invMass_y
diffRawMass = df.invMass_rawSC_x-df.invMass_rawSC_y

diffDf=df
diffDf['e1EnergyDiff']=diffE1
diffDf['e2EnergyDiff']=diffE2
diffDf['e1RawEnergyDiff']=diffRawE1
diffDf['e2RawEnergyDiff']=diffRawE2
diffDf['diffMass']=diffMass
diffDf['diffRawMass']=diffRawMass

#diffDf=diffDf.round(3)

#diffDf=diffDf[(abs(diffDf.e1RawEnergyDiff)>0.1) | (abs(diffDf.e2RawEnergyDiff)>0.1) | (abs(diffDf.diffRawMass)>0.1)]
diffDf=diffDf[diffDf.eventNumber==184638164]

#diffDf=diffDf.drop(columns=['e1Charge', 'e2Charge', 'e1IsDead_x', 'e1IsDead_y', 'e2IsDead_x', 'e2IsDead_y'])

#diffDf=diffDf[(diffDf.e1Eta_x!=diffDf.e1Eta_y) & (diffDf.e2Eta_x!=diffDf.e2Eta_y) & (diffDf.e1Phi_x!=diffDf.e1Phi_y)& (diffDf.e2Phi_x!=diffDf.e2Phi_y)]
#diffDf=diffDf[(diffDf.e1RawEnergy_x<diffDf.e1RawEnergy_y) | (diffDf.e2RawEnergy_x<diffDf.e2RawEnergy_y) ]

print list(diffDf.columns.values),  diffDf.shape[0]

#tmpDiffDf=diffDf[diffDf.diffRawMass!=0]
tmpDiffDf=diffDf
print tmpDiffDf[['runNumber', 'lumiBlock', 'eventNumber',#'e1Eta_x', 'e1Eta_y',  'e1Phi_x', 'e1Phi_y',  'e2Eta_x', 'e2Eta_y',  'e2Phi_x', 'e2Phi_y', 
#'e1EtaSC_x', 'e1EtaSC_y','e2EtaSC_x', 'e2EtaSC_y', 'e1PhiSC_x', 'e1PhiSC_y','e2PhiSC_x', 'e2PhiSC_y',# 'e1Energy_x', 
'e1RawEnergy_x', #'e2Energy_x', 
'e2RawEnergy_x', #'invMass_x', 
#                  'invMass_rawSC_x', 
   'e1RawEnergy_y', 'e2RawEnergy_y',# 'invMass_y', 
#                  'invMass_rawSC_y',# 'e1EnergyDiff', 'e2EnergyDiff',
#'e1IsRecovered', 'e2IsRecovered',
                  'e1RawEnergyDiff', 'e2RawEnergyDiff', 'diffRawMass']]#[["runNumber", "lumiBlock", "eventNumber", 'e1RawEnergyDiff', 'e2RawEnergyDiff', 'diffRawMass']]

#outfile=open("eventsToCheck.txt", "w")
#for x in range(0,len(diffDf)):
#    outfile.write( '%i:%i:%i\n' %(diffDf.iloc[x,:].runNumber, diffDf.iloc[x,:].lumiBlock, diffDf.iloc[x,:].eventNumber))

#outfile.close()

if createHistos:

    outrootfile=ROOT.TFile.Open("recoveryComparison_0_15.root", "RECREATE")
    outrootfile.cd()
    rInvMass=ROOT.TH1F("InvMass_15", "InvMass_15", 40, 0, 200)
    pInvMass=ROOT.TH1F("InvMass_0", "InvMass_0", 40, 0, 200)
    rInvMassRaw=ROOT.TH1F("InvMassRaw_15", "InvMassRaw_15", 40, 0, 200)
    pInvMassRaw=ROOT.TH1F("InvMassRaw_0", "InvMassRaw_0", 40, 0, 200)
    
    re1R9=ROOT.TH1F("e1R9_15", "e1R9_15", 24, 0, 1.2)
    pe1R9=ROOT.TH1F("e1R9_0", "e1R9_0", 24, 0, 1.2)
    re2R9=ROOT.TH1F("e2R9_15", "e2R9_15", 24, 0, 1.2)
    pe2R9=ROOT.TH1F("e2R9_0", "e2R9_0", 24, 0, 1.2)

    re1RawEnergy=ROOT.TH1F("e1RawEnergy_15", "e1RawEnergy_15", 40, 0, 200)
    pe1RawEnergy=ROOT.TH1F("e1RawEnergy_0", "e1RawEnergy_0", 40, 0, 200)
    re2RawEnergy=ROOT.TH1F("e2RawEnergy_15", "e2RawEnergy_15", 40, 0, 200)
    pe2RawEnergy=ROOT.TH1F("e2RawEnergy_0", "e2RawEnergy_0", 40, 0, 200)
    re1Energy=ROOT.TH1F("e1Energy_15", "e1Energy_15", 40, 0, 200)
    pe1Energy=ROOT.TH1F("e1Energy_0", "e1Energy_0", 40, 0, 200)
    re2Energy=ROOT.TH1F("e2Energy_15", "e2Energy_15", 40, 0, 200)
    pe2Energy=ROOT.TH1F("e2Energy_0", "e2Energy_0", 40, 0, 200)
    
    re1GenEnergy=ROOT.TH1F("e1GenEnergy_15", "e1GenEnergy_15", 40, 0, 200)
    pe1GenEnergy=ROOT.TH1F("e1GenEnergy_0", "e1GenEnergy_0", 40, 0, 200)
    re2GenEnergy=ROOT.TH1F("e2GenEnergy_15", "e2GenEnergy_15", 40, 0, 200)
    pe2GenEnergy=ROOT.TH1F("e2GenEnergy_0", "e2GenEnergy_0", 40, 0, 200)
    
    re1RawGenRatio=ROOT.TH1F("e1RawGenRatio_15", "e1RawGenRatio_15", 40, 0, 2)
    pe1RawGenRatio=ROOT.TH1F("e1RawGenRatio_0", "e1RawGenRatio_0", 40, 0, 2)
    re2RawGenRatio=ROOT.TH1F("e2RawGenRatio_15", "e2RawGenRatio_15", 40, 0, 2)
    pe2RawGenRatio=ROOT.TH1F("e2RawGenRatio_0", "e2RawGenRatio_0", 40, 0, 2)
    re1EnGenRatio=ROOT.TH1F("e1EnGenRatio_15", "e1EnGenRatio_15", 40, 0, 2)
    pe1EnGenRatio=ROOT.TH1F("e1EnGenRatio_0", "e1EnGenRatio_0", 40, 0, 2)
    re2EnGenRatio=ROOT.TH1F("e2EnGenRatio_15", "e2EnGenRatio_15", 40, 0, 2)
    pe2EnGenRatio=ROOT.TH1F("e2EnGenRatio_0", "e2EnGenRatio_0", 40, 0, 2)
    
    
    re1SigmaIetaIeta=ROOT.TH1F("e1SigmaIetaIeta_15","e1SigmaIetaIeta_15", 100, 0, 0.1)
    re2SigmaIetaIeta=ROOT.TH1F("e2SigmaIetaIeta_15","e2SigmaIetaIeta_15", 100, 0, 0.1)
    pe1SigmaIetaIeta=ROOT.TH1F("e1SigmaIetaIeta_0","e1SigmaIetaIeta_0", 100, 0, 0.1)
    pe2SigmaIetaIeta=ROOT.TH1F("e2SigmaIetaIeta_0","e2SigmaIetaIeta_0", 100, 0, 0.1)
    
    
    
    r_all_InvMass=ROOT.TH1F("InvMass_all_15", "InvMass_all_15", 40, 0, 200)
    p_all_InvMass=ROOT.TH1F("InvMass_all_0", "InvMass_all_0", 40, 0, 200)
    r_all_InvMassRaw=ROOT.TH1F("InvMassRaw_all_15", "InvMassRaw_all_15", 40, 0, 200)
    p_all_InvMassRaw=ROOT.TH1F("InvMassRaw_all_0", "InvMassRaw_all_0", 40, 0, 200)

    r_all_e1R9=ROOT.TH1F("e1R9_all_15", "e1R9_all_15", 24, 0, 1.2)
    p_all_e1R9=ROOT.TH1F("e1R9_all_0", "e1R9_all_0", 24, 0, 1.2)
    r_all_e2R9=ROOT.TH1F("e2R9_all_15", "e2R9_all_15", 24, 0, 1.2)
    p_all_e2R9=ROOT.TH1F("e2R9_all_0", "e2R9_all_0", 24, 0, 1.2)
    
    r_all_e1RawEnergy=ROOT.TH1F("e1RawEnergy_all_15", "e1RawEnergy_all_15", 40, 0, 200)
    p_all_e1RawEnergy=ROOT.TH1F("e1RawEnergy_all_0", "e1RawEnergy_all_0", 40, 0, 200)
    r_all_e2RawEnergy=ROOT.TH1F("e2RawEnergy_all_15", "e2RawEnergy_all_15", 40, 0, 200)
    p_all_e2RawEnergy=ROOT.TH1F("e2RawEnergy_all_0", "e2RawEnergy_all_0", 40, 0, 200)
    r_all_e1Energy=ROOT.TH1F("e1Energy_all_15", "e1Energy_all_15", 40, 0, 200)
    p_all_e1Energy=ROOT.TH1F("e1Energy_all_0", "e1Energy_all_0", 40, 0, 200)
    r_all_e2Energy=ROOT.TH1F("e2Energy_all_15", "e2Energy_all_15", 40, 0, 200)
    p_all_e2Energy=ROOT.TH1F("e2Energy_all_0", "e2Energy_all_0", 40, 0, 200)
    
    r_all_e1GenEnergy=ROOT.TH1F("e1GenEnergy_all_15", "e1GenEnergy_all_15", 40, 0, 200)
    p_all_e1GenEnergy=ROOT.TH1F("e1GenEnergy_all_0", "e1GenEnergy_all_0", 40, 0, 200)
    r_all_e2GenEnergy=ROOT.TH1F("e2GenEnergy_all_15", "e2GenEnergy_all_15", 40, 0, 200)
    p_all_e2GenEnergy=ROOT.TH1F("e2GenEnergy_all_0", "e2GenEnergy_all_0", 40, 0, 200)

    r_all_e1RawGenRatio=ROOT.TH1F("e1RawGenRatio_all_15", "e1RawGenRatio_all_15", 40, 0, 2)
    p_all_e1RawGenRatio=ROOT.TH1F("e1RawGenRatio_all_0", "e1RawGenRatio_all_0", 40, 0, 2)
    r_all_e2RawGenRatio=ROOT.TH1F("e2RawGenRatio_all_15", "e2RawGenRatio_all_15", 40, 0, 2)
    p_all_e2RawGenRatio=ROOT.TH1F("e2RawGenRatio_all_0", "e2RawGenRatio_all_0", 40, 0, 2)
    r_all_e1EnGenRatio=ROOT.TH1F("e1EnGenRatio_all_15", "e1EnGenRatio_all_15", 40, 0, 2)
    p_all_e1EnGenRatio=ROOT.TH1F("e1EnGenRatio_all_0", "e1EnGenRatio_all_0", 40, 0, 2)
    r_all_e2EnGenRatio=ROOT.TH1F("e2EnGenRatio_all_15", "e2EnGenRatio_all_15", 40, 0, 2)
    p_all_e2EnGenRatio=ROOT.TH1F("e2EnGenRatio_all_0", "e2EnGenRatio_all_0", 40, 0, 2)

    r_all_e1SigmaIetaIeta=ROOT.TH1F("e1SigmaIetaIeta_all_15","e1SigmaIetaIeta_all_15", 100, 0, 0.1)
    r_all_e2SigmaIetaIeta=ROOT.TH1F("e2SigmaIetaIeta_all_15","e2SigmaIetaIeta_all_15", 100, 0, 0.1)
    p_all_e1SigmaIetaIeta=ROOT.TH1F("e1SigmaIetaIeta_all_0","e1SigmaIetaIeta_all_0", 100, 0, 0.1)
    p_all_e2SigmaIetaIeta=ROOT.TH1F("e2SigmaIetaIeta_all_0","e2SigmaIetaIeta_all_0", 100, 0, 0.1)


    

    for x in range(0,len(diffDf)):

        rInvMass.Fill(diffDf.iloc[x,:].invMass_y)
        rInvMassRaw.Fill(diffDf.iloc[x,:].invMass_rawSC_y)
            
        pInvMass.Fill(diffDf.iloc[x,:].invMass_x)
        pInvMassRaw.Fill(diffDf.iloc[x,:].invMass_rawSC_x)
        if (diffDf.iloc[x,:].e1IsRecovered):
            re1R9.Fill(diffDf.iloc[x,:].e1R9_y)
            re1RawEnergy.Fill(diffDf.iloc[x,:].e1RawEnergy_y)
            re1Energy.Fill(diffDf.iloc[x,:].e1Energy_y)
            ##re1GenEnergy.Fill(diffDf.iloc[x,:].e1GenEnergy)
            ##re1RawGenRatio.Fill(diffDf.iloc[x,:].e1RawEnergy_y/diffDf.iloc[x,:].e1GenEnergy)
            ##re1EnGenRatio.Fill(diffDf.iloc[x,:].e1Energy_y/diffDf.iloc[x,:].e1GenEnergy)
            re1SigmaIetaIeta.Fill(diffDf.iloc[x,:].e1SigmaIetaIeta_y)
            pe1R9.Fill(diffDf.iloc[x,:].e1R9_x)
            pe1RawEnergy.Fill(diffDf.iloc[x,:].e1RawEnergy_x)
            pe1Energy.Fill(diffDf.iloc[x,:].e1Energy_x)
            ##pe1GenEnergy.Fill(diffDf.iloc[x,:].e1GenEnergy)
            ##pe1RawGenRatio.Fill(diffDf.iloc[x,:].e1RawEnergy_x/diffDf.iloc[x,:].e1GenEnergy)
            ##pe1EnGenRatio.Fill(diffDf.iloc[x,:].e1Energy_x/diffDf.iloc[x,:].e1GenEnergy)
            ##pe1SigmaIetaIeta.Fill(diffDf.iloc[x,:].e1SigmaIetaIeta_x)
            
        elif(diffDf.iloc[x,:].e2IsRecovered):
            re2R9.Fill(diffDf.iloc[x,:].e2R9_y)
            re2RawEnergy.Fill(diffDf.iloc[x,:].e2RawEnergy_y)
            re2Energy.Fill(diffDf.iloc[x,:].e2Energy_y)
            ##re2GenEnergy.Fill(diffDf.iloc[x,:].e2GenEnergy)
            ##re2RawGenRatio.Fill(diffDf.iloc[x,:].e2RawEnergy_y/diffDf.iloc[x,:].e2GenEnergy)
            ##re2EnGenRatio.Fill(diffDf.iloc[x,:].e2Energy_y/diffDf.iloc[x,:].e2GenEnergy)
            re2SigmaIetaIeta.Fill(diffDf.iloc[x,:].e2SigmaIetaIeta_y)
            pe2R9.Fill(diffDf.iloc[x,:].e2R9_x)
            pe2RawEnergy.Fill(diffDf.iloc[x,:].e2RawEnergy_x)
            pe2Energy.Fill(diffDf.iloc[x,:].e2Energy_x)
            ##pe2GenEnergy.Fill(diffDf.iloc[x,:].e2GenEnergy)
            ##pe2RawGenRatio.Fill(diffDf.iloc[x,:].e2RawEnergy_x/diffDf.iloc[x,:].e2GenEnergy)
            ##pe2EnGenRatio.Fill(diffDf.iloc[x,:].e2Energy_x/diffDf.iloc[x,:].e2GenEnergy)
            pe2SigmaIetaIeta.Fill(diffDf.iloc[x,:].e2SigmaIetaIeta_x)

    for x in range(0,len(indata0)):
        p_all_InvMass.Fill(indata0.iloc[x,:].invMass)
        p_all_InvMassRaw.Fill(indata0.iloc[x,:].invMass_rawSC)
        if (indata0.iloc[x,:].e1IsRecovered):
            p_all_e1R9.Fill(indata0.iloc[x,:].e1R9)
            p_all_e1RawEnergy.Fill(indata0.iloc[x,:].e1RawEnergy)
            p_all_e1Energy.Fill(indata0.iloc[x,:].e1Energy)
            ##p_all_e1GenEnergy.Fill(indata0.iloc[x,:].e1GenEnergy)
            ##p_all_e1RawGenRatio.Fill(indata0.iloc[x,:].e1RawEnergy/indata0.iloc[x,:].e1GenEnergy)
            ##p_all_e1EnGenRatio.Fill(indata0.iloc[x,:].e1Energy/indata0.iloc[x,:].e1GenEnergy)
            p_all_e1SigmaIetaIeta.Fill(indata0.iloc[x,:].e1SigmaIetaIeta)

        elif(indata0.iloc[x,:].e2IsRecovered):
            p_all_e2R9.Fill(indata0.iloc[x,:].e2R9)
            p_all_e2RawEnergy.Fill(indata0.iloc[x,:].e2RawEnergy)
            p_all_e2Energy.Fill(indata0.iloc[x,:].e2Energy)
            ##p_all_e2GenEnergy.Fill(indata0.iloc[x,:].e2GenEnergy)
            ##p_all_e2RawGenRatio.Fill(indata0.iloc[x,:].e2RawEnergy/indata0.iloc[x,:].e2GenEnergy)
            ##p_all_e2EnGenRatio.Fill(indata0.iloc[x,:].e2Energy/indata0.iloc[x,:].e2GenEnergy)
            p_all_e2SigmaIetaIeta.Fill(indata0.iloc[x,:].e2SigmaIetaIeta)

    for x in range(0,len(indata_noCh)):
        r_all_InvMass.Fill(indata_noCh.iloc[x,:].invMass)
        r_all_InvMassRaw.Fill(indata_noCh.iloc[x,:].invMass_rawSC)
        if (indata_noCh.iloc[x,:].e1IsRecovered):
            r_all_e1R9.Fill(indata_noCh.iloc[x,:].e1R9)
            r_all_e1RawEnergy.Fill(indata_noCh.iloc[x,:].e1RawEnergy)
            r_all_e1Energy.Fill(indata_noCh.iloc[x,:].e1Energy)
            ##r_all_e1GenEnergy.Fill(indata_noCh.iloc[x,:].e1GenEnergy)
            ##r_all_e1RawGenRatio.Fill(indata_noCh.iloc[x,:].e1RawEnergy/indata_noCh.iloc[x,:].e1GenEnergy)
            ##r_all_e1EnGenRatio.Fill(indata_noCh.iloc[x,:].e1Energy/indata_noCh.iloc[x,:].e1GenEnergy)
            r_all_e1SigmaIetaIeta.Fill(indata_noCh.iloc[x,:].e1SigmaIetaIeta)

        elif(indata_noCh.iloc[x,:].e2IsRecovered):
            r_all_e2R9.Fill(indata_noCh.iloc[x,:].e2R9)
            r_all_e2RawEnergy.Fill(indata_noCh.iloc[x,:].e2RawEnergy)
            r_all_e2Energy.Fill(indata_noCh.iloc[x,:].e2Energy)
            ##r_all_e2GenEnergy.Fill(indata_noCh.iloc[x,:].e2GenEnergy)
            ##r_all_e2RawGenRatio.Fill(indata_noCh.iloc[x,:].e2RawEnergy/indata_noCh.iloc[x,:].e2GenEnergy)
            ##r_all_e2EnGenRatio.Fill(indata_noCh.iloc[x,:].e2Energy/indata_noCh.iloc[x,:].e2GenEnergy)
            r_all_e2SigmaIetaIeta.Fill(indata_noCh.iloc[x,:].e2SigmaIetaIeta)
            

    

    outrootfile.cd()
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

    p_all_InvMass.Write()
    p_all_InvMassRaw.Write()
    p_all_e1R9.Write()
    p_all_e2R9.Write()
    p_all_e1RawEnergy.Write()
    p_all_e2RawEnergy.Write()
    p_all_e1Energy.Write()
    p_all_e2Energy.Write()
    p_all_e1SigmaIetaIeta.Write()
    p_all_e2SigmaIetaIeta.Write()
    p_all_e2GenEnergy.Write()
    p_all_e2RawGenRatio.Write()
    p_all_e2EnGenRatio.Write()
    p_all_e1GenEnergy.Write()
    p_all_e1RawGenRatio.Write()
    p_all_e1EnGenRatio.Write()

    r_all_InvMass.Write()
    r_all_InvMassRaw.Write()
    r_all_e1R9.Write()
    r_all_e2R9.Write()
    r_all_e1RawEnergy.Write()
    r_all_e2RawEnergy.Write()
    r_all_e1Energy.Write()
    r_all_e2Energy.Write()
    r_all_e1SigmaIetaIeta.Write()
    r_all_e2SigmaIetaIeta.Write()
    r_all_e2GenEnergy.Write()
    r_all_e2RawGenRatio.Write()
    r_all_e2EnGenRatio.Write()
    r_all_e1GenEnergy.Write()
    r_all_e1RawGenRatio.Write()
    r_all_e1EnGenRatio.Write()
    
    outrootfile.Close()

