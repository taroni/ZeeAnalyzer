import ROOT
import uproot
import numpy as np
import pandas as pd
import array as arr

import argparse


def deltaR(eta1, phi1, eta2, phi2):
    deltaEta=abs(eta1-eta2)
    deltaPhi=abs(phi1-phi2) 
    #   print deltaPhi
    deltaPhi[deltaPhi>np.pi]=deltaPhi-2*np.pi
    #print deltaPhi
    return np.sqrt(deltaEta**2+deltaPhi**2)


parser = argparse.ArgumentParser(description='createHistograms')
parser.add_argument('--createHistograms', dest='histos', action='store', default=False, help='boolean if saving histos', required=False)

args = parser.parse_args()
createHistos=bool(args.histos)

tree0 =uproot.open("electronTree_PierreTag_MCPU_0GeV_reweight.root")['ntupler']['selected']
tree_noCh=uproot.open("electronTree_PierreTag_MCPU_0GeV_noRecov_reweight.root")['ntupler']['selected']

indata0=tree0.pandas.df(['runNumber', 'lumiBlock', 'eventNumber','e1Charge','e2Charge','e1Eta', 'e1Phi', 'e1R9', 'e2Eta', 'e2Phi', 'e2R9', 'e1SeedIEta', 'e1SeedIPhi', 'e2SeedIEta', 'e2SeedIPhi', 'e1SigmaIetaIeta', 'e2SigmaIetaIeta', 'e1EtaSC', 'e2EtaSC', 'e1PhiSC', 'e2PhiSC',   'e1Energy', 'e1RawEnergy', 'e2Energy', 'e2RawEnergy',  'invMass', 'invMass_rawSC', 'e1IsDead', 'e1IsRecovered', 'e2IsDead', 'e2IsRecovered'
])
indata_noCh=tree_noCh.pandas.df(['runNumber', 'lumiBlock', 'eventNumber', 'e1Charge', 'e2Charge', 'e1Eta', 'e1Phi', 'e1R9', 'e2Eta', 'e2Phi', 'e2R9', 'e1SeedIEta', 'e1SeedIPhi', 'e2SeedIEta', 'e2SeedIPhi', 'e1SigmaIetaIeta', 'e2SigmaIetaIeta', 'e1EtaSC', 'e2EtaSC', 'e1PhiSC', 'e2PhiSC', 'e1Energy', 'e1RawEnergy', 'e2Energy', 'e2RawEnergy', 'invMass', 'invMass_rawSC', 'e1IsDead', 'e1IsRecovered', 'e2IsDead', 'e2IsRecovered'
])


indata_noCh=indata_noCh[indata_noCh.e1Charge*indata_noCh.e2Charge==-1]
indata0=indata0[indata0.e1Charge*indata0.e2Charge==-1]

indata0=indata0[(indata0.e1IsRecovered==1) | (indata0.e2IsRecovered==1)]

evtlist0 = set(["%i:%i:%i" %(indata0.iloc[x,:].runNumber, indata0.iloc[x,:].lumiBlock, indata0.iloc[x,:].eventNumber) for x in range(0, len(indata0))])
evtlist = set(["%i:%i:%i" %(indata_noCh.iloc[x,:].runNumber, indata_noCh.iloc[x,:].lumiBlock, indata_noCh.iloc[x,:].eventNumber) for x in range(0, len(indata_noCh))])
print indata_noCh.head(20)
print indata0.head(20)

#keep only events in both datasets
df=pd.merge(indata_noCh, indata0, how='inner',on=["runNumber", "lumiBlock", "eventNumber"], indicator=True) 

#keep only events with matching electron seeds
#df=df[(abs(df.e1SeedIPhi_x-df.e1SeedIPhi_y)<2)&(abs(df.e1SeedIEta_x-df.e1SeedIEta_y)<2) &(abs(df.e2SeedIPhi_x-df.e2SeedIPhi_y)<2)&(abs(df.e2SeedIEta_x-df.e2SeedIEta_y)<2)]
print list(df.columns.values)

df=df[(abs(df.e1SeedIPhi_x-df.e1SeedIPhi_y)<5)&(abs(df.e1SeedIEta_x-df.e1SeedIEta_y)<5) &(abs(df.e2SeedIPhi_x-df.e2SeedIPhi_y)<5)&(abs(df.e2SeedIEta_x-df.e2SeedIEta_y)<5) | (abs(df.e1SeedIPhi_x-df.e2SeedIPhi_y)<5)&(abs(df.e1SeedIEta_x-df.e2SeedIEta_y)<5) &(abs(df.e2SeedIPhi_x-df.e1SeedIPhi_y)<5)&(abs(df.e2SeedIEta_x-df.e1SeedIEta_y)<5)]
print df.shape[0]
df=df[((deltaR(df.e1Eta_x, df.e1Phi_x, df.e1Eta_y, df.e1Phi_y)<0.5) & (deltaR(df.e2Eta_x, df.e2Phi_x, df.e2Eta_y, df.e2Phi_y)<0.5)) | ((deltaR(df.e1Eta_x, df.e1Phi_x, df.e2Eta_y, df.e2Phi_y)<0.5) & (deltaR(df.e2Eta_x, df.e2Phi_x, df.e1Eta_y, df.e1Phi_y)<0.5))]
#print df
#print list(df.columns.values)
print df.loc[:, df.columns.isin([ 'e1SeedIEta_x', 'e1SeedIPhi_x', 'e2SeedIEta_x', 'e2SeedIPhi_x','e1SeedIEta_y', 'e1SeedIPhi_y', 'e2SeedIEta_y', 'e2SeedIPhi_y', 'invMass_rawSC_x','invMass_rawSC_y'])]
print df.loc[:, df.columns.isin([ 'e1Eta_x', 'e1Phi_x', 'e2Eta_x', 'e2Phi_x','e1Eta_y', 'e1Phi_y', 'e2Eta_y', 'e2Phi_y', 'invMass_rawSC_x','invMass_rawSC_y'])]

print df.shape[0]

if createHistos:

    outrootfile=ROOT.TFile.Open("recoveryComparisonMC_noRec_rec.root", "RECREATE")
    outrootfile.cd()
    rInvMass=ROOT.TH1F("InvMass_rec", "InvMass_rec", 40, 0, 200)
    pInvMass=ROOT.TH1F("InvMass_noRec", "InvMass_noRec", 40, 0, 200)
    rInvMassRaw=ROOT.TH1F("InvMassRaw_rec", "InvMassRaw_rec", 40, 0, 200)
    pInvMassRaw=ROOT.TH1F("InvMassRaw_noRec", "InvMassRaw_noRec", 40, 0, 200)
    
    re1R9=ROOT.TH1F("e1R9_rec", "e1R9_rec", 24, 0, 1.2)
    pe1R9=ROOT.TH1F("e1R9_noRec", "e1R9_noRec", 24, 0, 1.2)
    re2R9=ROOT.TH1F("e2R9_rec", "e2R9_rec", 24, 0, 1.2)
    pe2R9=ROOT.TH1F("e2R9_noRec", "e2R9_noRec", 24, 0, 1.2)

    re1RawEnergy=ROOT.TH1F("e1RawEnergy_rec", "e1RawEnergy_rec", 40, 0, 200)
    pe1RawEnergy=ROOT.TH1F("e1RawEnergy_noRec", "e1RawEnergy_noRec", 40, 0, 200)
    re2RawEnergy=ROOT.TH1F("e2RawEnergy_rec", "e2RawEnergy_rec", 40, 0, 200)
    pe2RawEnergy=ROOT.TH1F("e2RawEnergy_noRec", "e2RawEnergy_noRec", 40, 0, 200)
    re1Energy=ROOT.TH1F("e1Energy_rec", "e1Energy_rec", 40, 0, 200)
    pe1Energy=ROOT.TH1F("e1Energy_noRec", "e1Energy_noRec", 40, 0, 200)
    re2Energy=ROOT.TH1F("e2Energy_rec", "e2Energy_rec", 40, 0, 200)
    pe2Energy=ROOT.TH1F("e2Energy_noRec", "e2Energy_noRec", 40, 0, 200)
    
    re1GenEnergy=ROOT.TH1F("e1GenEnergy_rec", "e1GenEnergy_rec", 40, 0, 200)
    pe1GenEnergy=ROOT.TH1F("e1GenEnergy_noRec", "e1GenEnergy_noRec", 40, 0, 200)
    re2GenEnergy=ROOT.TH1F("e2GenEnergy_rec", "e2GenEnergy_rec", 40, 0, 200)
    pe2GenEnergy=ROOT.TH1F("e2GenEnergy_noRec", "e2GenEnergy_noRec", 40, 0, 200)
    
    re1RawGenRatio=ROOT.TH1F("e1RawGenRatio_rec", "e1RawGenRatio_rec", 40, 0, 2)
    pe1RawGenRatio=ROOT.TH1F("e1RawGenRatio_noRec", "e1RawGenRatio_noRec", 40, 0, 2)
    re2RawGenRatio=ROOT.TH1F("e2RawGenRatio_rec", "e2RawGenRatio_rec", 40, 0, 2)
    pe2RawGenRatio=ROOT.TH1F("e2RawGenRatio_noRec", "e2RawGenRatio_noRec", 40, 0, 2)
    re1EnGenRatio=ROOT.TH1F("e1EnGenRatio_rec", "e1EnGenRatio_rec", 40, 0, 2)
    pe1EnGenRatio=ROOT.TH1F("e1EnGenRatio_noRec", "e1EnGenRatio_noRec", 40, 0, 2)
    re2EnGenRatio=ROOT.TH1F("e2EnGenRatio_rec", "e2EnGenRatio_rec", 40, 0, 2)
    pe2EnGenRatio=ROOT.TH1F("e2EnGenRatio_noRec", "e2EnGenRatio_noRec", 40, 0, 2)
    
    
    re1SigmaIetaIeta=ROOT.TH1F("e1SigmaIetaIeta_rec","e1SigmaIetaIeta_rec", 100, 0, 0.1)
    re2SigmaIetaIeta=ROOT.TH1F("e2SigmaIetaIeta_rec","e2SigmaIetaIeta_rec", 100, 0, 0.1)
    pe1SigmaIetaIeta=ROOT.TH1F("e1SigmaIetaIeta_noRec","e1SigmaIetaIeta_noRec", 100, 0, 0.1)
    pe2SigmaIetaIeta=ROOT.TH1F("e2SigmaIetaIeta_noRec","e2SigmaIetaIeta_noRec", 100, 0, 0.1)
    
    
    
    r_all_InvMass=ROOT.TH1F("InvMass_all_rec", "InvMass_all_rec", 40, 0, 200)
    p_all_InvMass=ROOT.TH1F("InvMass_all_noRec", "InvMass_all_noRec", 40, 0, 200)
    r_all_InvMassRaw=ROOT.TH1F("InvMassRaw_all_rec", "InvMassRaw_all_rec", 40, 0, 200)
    p_all_InvMassRaw=ROOT.TH1F("InvMassRaw_all_noRec", "InvMassRaw_all_noRec", 40, 0, 200)

    r_all_e1R9=ROOT.TH1F("e1R9_all_rec", "e1R9_all_rec", 24, 0, 1.2)
    p_all_e1R9=ROOT.TH1F("e1R9_all_noRec", "e1R9_all_noRec", 24, 0, 1.2)
    r_all_e2R9=ROOT.TH1F("e2R9_all_rec", "e2R9_all_rec", 24, 0, 1.2)
    p_all_e2R9=ROOT.TH1F("e2R9_all_noRec", "e2R9_all_noRec", 24, 0, 1.2)
    
    r_all_e1RawEnergy=ROOT.TH1F("e1RawEnergy_all_rec", "e1RawEnergy_all_rec", 40, 0, 200)
    p_all_e1RawEnergy=ROOT.TH1F("e1RawEnergy_all_noRec", "e1RawEnergy_all_noRec", 40, 0, 200)
    r_all_e2RawEnergy=ROOT.TH1F("e2RawEnergy_all_rec", "e2RawEnergy_all_rec", 40, 0, 200)
    p_all_e2RawEnergy=ROOT.TH1F("e2RawEnergy_all_noRec", "e2RawEnergy_all_noRec", 40, 0, 200)
    r_all_e1Energy=ROOT.TH1F("e1Energy_all_rec", "e1Energy_all_rec", 40, 0, 200)
    p_all_e1Energy=ROOT.TH1F("e1Energy_all_noRec", "e1Energy_all_noRec", 40, 0, 200)
    r_all_e2Energy=ROOT.TH1F("e2Energy_all_rec", "e2Energy_all_rec", 40, 0, 200)
    p_all_e2Energy=ROOT.TH1F("e2Energy_all_noRec", "e2Energy_all_noRec", 40, 0, 200)
    
    r_all_e1GenEnergy=ROOT.TH1F("e1GenEnergy_all_rec", "e1GenEnergy_all_rec", 40, 0, 200)
    p_all_e1GenEnergy=ROOT.TH1F("e1GenEnergy_all_noRec", "e1GenEnergy_all_noRec", 40, 0, 200)
    r_all_e2GenEnergy=ROOT.TH1F("e2GenEnergy_all_rec", "e2GenEnergy_all_rec", 40, 0, 200)
    p_all_e2GenEnergy=ROOT.TH1F("e2GenEnergy_all_noRec", "e2GenEnergy_all_noRec", 40, 0, 200)

    r_all_e1RawGenRatio=ROOT.TH1F("e1RawGenRatio_all_rec", "e1RawGenRatio_all_rec", 40, 0, 2)
    p_all_e1RawGenRatio=ROOT.TH1F("e1RawGenRatio_all_noRec", "e1RawGenRatio_all_noRec", 40, 0, 2)
    r_all_e2RawGenRatio=ROOT.TH1F("e2RawGenRatio_all_rec", "e2RawGenRatio_all_rec", 40, 0, 2)
    p_all_e2RawGenRatio=ROOT.TH1F("e2RawGenRatio_all_noRec", "e2RawGenRatio_all_noRec", 40, 0, 2)
    r_all_e1EnGenRatio=ROOT.TH1F("e1EnGenRatio_all_rec", "e1EnGenRatio_all_rec", 40, 0, 2)
    p_all_e1EnGenRatio=ROOT.TH1F("e1EnGenRatio_all_noRec", "e1EnGenRatio_all_noRec", 40, 0, 2)
    r_all_e2EnGenRatio=ROOT.TH1F("e2EnGenRatio_all_rec", "e2EnGenRatio_all_rec", 40, 0, 2)
    p_all_e2EnGenRatio=ROOT.TH1F("e2EnGenRatio_all_noRec", "e2EnGenRatio_all_noRec", 40, 0, 2)

    r_all_e1SigmaIetaIeta=ROOT.TH1F("e1SigmaIetaIeta_all_rec","e1SigmaIetaIeta_all_rec", 100, 0, 0.1)
    r_all_e2SigmaIetaIeta=ROOT.TH1F("e2SigmaIetaIeta_all_rec","e2SigmaIetaIeta_all_rec", 100, 0, 0.1)
    p_all_e1SigmaIetaIeta=ROOT.TH1F("e1SigmaIetaIeta_all_noRec","e1SigmaIetaIeta_all_noRec", 100, 0, 0.1)
    p_all_e2SigmaIetaIeta=ROOT.TH1F("e2SigmaIetaIeta_all_noRec","e2SigmaIetaIeta_all_noRec", 100, 0, 0.1)


    r_allrecov_e1R9=ROOT.TH1F("e1R9_allrecov_rec", "e1R9_allrecov_rec", 24, 0, 1.2)
    r_allrecov_e2R9=ROOT.TH1F("e2R9_allrecov_rec", "e2R9_allrecov_rec", 24, 0, 1.2)
    r_allrecov_e1RawEnergy=ROOT.TH1F("e1RawEnergy_allrecov_rec", "e1RawEnergy_allrecov_rec", 40, 0, 200)
    r_allrecov_e2RawEnergy=ROOT.TH1F("e2RawEnergy_allrecov_rec", "e2RawEnergy_allrecov_rec", 40, 0, 200)
    r_allrecov_e1Energy=ROOT.TH1F("e1Energy_allrecov_rec", "e1Energy_allrecov_rec", 40, 0, 200)
    r_allrecov_e2Energy=ROOT.TH1F("e2Energy_allrecov_rec", "e2Energy_allrecov_rec", 40, 0, 200)
    r_allrecov_e1SigmaIetaIeta=ROOT.TH1F("e1SigmaIetaIeta_allrecov_rec","e1SigmaIetaIeta_allrecov_rec", 100, 0, 0.1)
    r_allrecov_e2SigmaIetaIeta=ROOT.TH1F("e2SigmaIetaIeta_allrecov_rec","e2SigmaIetaIeta_allrecov_rec", 100, 0, 0.1)


    ## entries ending with _y are the variables from recovered databased

    for x in range(0,len(df)):

        rInvMass.Fill(df.iloc[x,:].invMass_y)
        rInvMassRaw.Fill(df.iloc[x,:].invMass_rawSC_y)
            
        pInvMass.Fill(df.iloc[x,:].invMass_x)
        pInvMassRaw.Fill(df.iloc[x,:].invMass_rawSC_x)
        if (df.iloc[x,:].e1IsRecovered_y):
            re1R9.Fill(df.iloc[x,:].e1R9_y)
            re1RawEnergy.Fill(df.iloc[x,:].e1RawEnergy_y)
            re1Energy.Fill(df.iloc[x,:].e1Energy_y)
            ##re1GenEnergy.Fill(df.iloc[x,:].e1GenEnergy)
            ##re1RawGenRatio.Fill(df.iloc[x,:].e1RawEnergy_y/df.iloc[x,:].e1GenEnergy)
            ##re1EnGenRatio.Fill(df.iloc[x,:].e1Energy_y/df.iloc[x,:].e1GenEnergy)
            re1SigmaIetaIeta.Fill(df.iloc[x,:].e1SigmaIetaIeta_y)
            pe1R9.Fill(df.iloc[x,:].e1R9_x)
            pe1RawEnergy.Fill(df.iloc[x,:].e1RawEnergy_x)
            pe1Energy.Fill(df.iloc[x,:].e1Energy_x)
            ##pe1GenEnergy.Fill(df.iloc[x,:].e1GenEnergy)
            ##pe1RawGenRatio.Fill(df.iloc[x,:].e1RawEnergy_x/df.iloc[x,:].e1GenEnergy)
            ##pe1EnGenRatio.Fill(df.iloc[x,:].e1Energy_x/df.iloc[x,:].e1GenEnergy)
            ##pe1SigmaIetaIeta.Fill(df.iloc[x,:].e1SigmaIetaIeta_x)
            
        elif(df.iloc[x,:].e2IsRecovered_y):
            re2R9.Fill(df.iloc[x,:].e2R9_y)
            re2RawEnergy.Fill(df.iloc[x,:].e2RawEnergy_y)
            re2Energy.Fill(df.iloc[x,:].e2Energy_y)
            ##re2GenEnergy.Fill(df.iloc[x,:].e2GenEnergy)
            ##re2RawGenRatio.Fill(df.iloc[x,:].e2RawEnergy_y/df.iloc[x,:].e2GenEnergy)
            ##re2EnGenRatio.Fill(df.iloc[x,:].e2Energy_y/df.iloc[x,:].e2GenEnergy)
            re2SigmaIetaIeta.Fill(df.iloc[x,:].e2SigmaIetaIeta_y)
            pe2R9.Fill(df.iloc[x,:].e2R9_x)
            pe2RawEnergy.Fill(df.iloc[x,:].e2RawEnergy_x)
            pe2Energy.Fill(df.iloc[x,:].e2Energy_x)
            ##pe2GenEnergy.Fill(df.iloc[x,:].e2GenEnergy)
            ##pe2RawGenRatio.Fill(df.iloc[x,:].e2RawEnergy_x/df.iloc[x,:].e2GenEnergy)
            ##pe2EnGenRatio.Fill(df.iloc[x,:].e2Energy_x/df.iloc[x,:].e2GenEnergy)
            pe2SigmaIetaIeta.Fill(df.iloc[x,:].e2SigmaIetaIeta_x)

    for x in range(0,len(indata_noCh)):
        p_all_InvMass.Fill(indata_noCh.iloc[x,:].invMass)
        p_all_InvMassRaw.Fill(indata_noCh.iloc[x,:].invMass_rawSC)
        p_all_e1R9.Fill(indata_noCh.iloc[x,:].e1R9)
        p_all_e1RawEnergy.Fill(indata_noCh.iloc[x,:].e1RawEnergy)
        p_all_e1Energy.Fill(indata_noCh.iloc[x,:].e1Energy)
        ##p_all_e1GenEnergy.Fill(indata_noCh.iloc[x,:].e1GenEnergy)
        ##p_all_e1RawGenRatio.Fill(indata_noCh.iloc[x,:].e1RawEnergy/indata_noCh.iloc[x,:].e1GenEnergy)
        ##p_all_e1EnGenRatio.Fill(indata_noCh.iloc[x,:].e1Energy/indata_noCh.iloc[x,:].e1GenEnergy)
        p_all_e1SigmaIetaIeta.Fill(indata_noCh.iloc[x,:].e1SigmaIetaIeta)

        p_all_e2R9.Fill(indata_noCh.iloc[x,:].e2R9)
        p_all_e2RawEnergy.Fill(indata_noCh.iloc[x,:].e2RawEnergy)
        p_all_e2Energy.Fill(indata_noCh.iloc[x,:].e2Energy)
        ##p_all_e2GenEnergy.Fill(indata_noCh.iloc[x,:].e2GenEnergy)
        ##p_all_e2RawGenRatio.Fill(indata_noCh.iloc[x,:].e2RawEnergy/indata_noCh.iloc[x,:].e2GenEnergy)
        ##p_all_e2EnGenRatio.Fill(indata_noCh.iloc[x,:].e2Energy/indata_noCh.iloc[x,:].e2GenEnergy)
        p_all_e2SigmaIetaIeta.Fill(indata_noCh.iloc[x,:].e2SigmaIetaIeta)

    for x in range(0,len(indata0)):
        r_all_InvMass.Fill(indata0.iloc[x,:].invMass)
        r_all_InvMassRaw.Fill(indata0.iloc[x,:].invMass_rawSC)
        r_all_e1R9.Fill(indata0.iloc[x,:].e1R9)
        r_all_e1RawEnergy.Fill(indata0.iloc[x,:].e1RawEnergy)
        r_all_e1Energy.Fill(indata0.iloc[x,:].e1Energy)
        r_all_e1SigmaIetaIeta.Fill(indata0.iloc[x,:].e1SigmaIetaIeta)
        r_all_e2R9.Fill(indata0.iloc[x,:].e2R9)
        r_all_e2RawEnergy.Fill(indata0.iloc[x,:].e2RawEnergy)
        r_all_e2Energy.Fill(indata0.iloc[x,:].e2Energy)
        r_all_e2SigmaIetaIeta.Fill(indata0.iloc[x,:].e2SigmaIetaIeta)
 
        if (indata0.iloc[x,:].e1IsRecovered):
            r_allrecov_e1R9.Fill(indata0.iloc[x,:].e1R9)
            r_allrecov_e1RawEnergy.Fill(indata0.iloc[x,:].e1RawEnergy)
            r_allrecov_e1Energy.Fill(indata0.iloc[x,:].e1Energy)
            r_allrecov_e1SigmaIetaIeta.Fill(indata0.iloc[x,:].e1SigmaIetaIeta)

        elif(indata0.iloc[x,:].e2IsRecovered):
            r_allrecov_e2R9.Fill(indata0.iloc[x,:].e2R9)
            r_allrecov_e2RawEnergy.Fill(indata0.iloc[x,:].e2RawEnergy)
            r_allrecov_e2Energy.Fill(indata0.iloc[x,:].e2Energy)
            r_allrecov_e2SigmaIetaIeta.Fill(indata0.iloc[x,:].e2SigmaIetaIeta)
            

    

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

    r_allrecov_e1R9.Write()
    r_allrecov_e2R9.Write()
    r_allrecov_e1RawEnergy.Write()
    r_allrecov_e2RawEnergy.Write()
    r_allrecov_e1Energy.Write()
    r_allrecov_e2Energy.Write()
    r_allrecov_e1SigmaIetaIeta.Write()
    r_allrecov_e2SigmaIetaIeta.Write()
    
    outrootfile.Close()

