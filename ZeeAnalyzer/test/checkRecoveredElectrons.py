import ROOT
import uproot
import numpy as np
import pandas as pd

tree0 =uproot.open("electronTreeZEE_PierreTag_sum8gt20_xtalInfo.root")['ntupler']['selected']

indata0=tree0.pandas.df(['runNumber', 'lumiBlock', 'eventNumber', 'e1Charge', 'e1Eta', 'e1Phi', 'e1R9', 'e2Charge', 'e2Eta', 'e2Phi', 'e2R9', 'e1SigmaIetaIeta', 'e2SigmaIetaIeta', 'e1GenEnergy', 'e2GenEnergy', 'e1EtaSC', 'e1PhiSC', 'e1Energy', 'e1RawEnergy', 'e2EtaSC', 'e2PhiSC', 'e2Energy', 'e2RawEnergy', 'invMass', 'invMass_rawSC', 'e1IsDead', 'e1IsRecovered', 'e2IsDead', 'e2IsRecovered', 'e1SeedIEta', 'e2SeedIEta', 'e1SeedIPhi', 'e2SeedIPhi'])

indata0=indata0[(indata0.e1IsRecovered==True) | ( indata0.e2IsRecovered==True)]
print indata0.shape[0]

indata0=indata0[(indata0.e1IsRecovered==True) & ( indata0.e2IsRecovered==True)]

print indata0.shape[0]
