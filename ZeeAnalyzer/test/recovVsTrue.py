import ROOT
import pandas as pd
import uproot
import numpy as np

def zero_pad(jagged, pad_size):
    ret = np.zeros((len(jagged), pad_size))
    for idx in xrange(len(jagged)):
        limit = min(len(jagged[idx]),pad_size)
        #print idx, limit, jagged[idx], ret[idx]
        ret[idx,:limit] = jagged[idx]
        #print ret[idx]
    return ret

tree0 =uproot.open("electronTreeZEE_PierreTag_sum8gt20_xtalInfo.root")['ntupler']['selected']
tree_noCh=uproot.open("electronTreeZEE_NoChange_sum8gt20_xtalInfo.root")['ntupler']['selected']

#['runNumber', 'lumiBlock', 'eventNumber', 'eventTime', 'nBX', 'nTruePU', 'puweight', 'e1Charge', 'e1Eta', 'e1Phi', 'e1R9', 'e2Charge', 'e2Eta', 'e2Phi', 'e2R9', 'e1SigmaIetaIeta', 'e2SigmaIetaIeta', 'e1GenEnergy', 'e2GenEnergy', 'e1EtaSC', 'e1PhiSC', 'e1Energy', 'e1RawEnergy', 'e2EtaSC', 'e2PhiSC', 'e2Energy', 'e2RawEnergy', 'invMass', 'invMass_rawSC', 'e1IsDead', 'e1IsRecovered', 'e2IsDead', 'e2IsRecovered', 'e1SeedIEta', 'e2SeedIEta', 'e1SeedIPhi', 'e2SeedIPhi', 'kGood', 'kPoorReco', 'kOutOfTimE', 'kFaultyHardware', 'kNoisy', 'kPoorCalib', 'kSaturated', 'kLeadingEdgeRecovered', 'kNeighboursRecovered', 'kTowerRecovered', 'kDead', 'kKilled', 'kTPSaturated', 'kL1SpikeFlag', 'kWeird', 'kDiWeird', 'kHasSwitchToGain6', 'kHasSwitchToGain1', 'kUnknown', 'xtalIeta', 'xtalIphi', 'xtalIsm', 'xtalIc', 'xtalRawId', 'xtalEn', 'vSum8', 'vAve', 'vNeigh0', 'vNeigh1', 'vNeigh2', 'vNeigh3', 'vNeigh5', 'vNeigh6', 'vNeigh7', 'vNeigh8']

myTree0=tree0.arrays(['runNumber', 'lumiBlock', 'eventNumber', 'xtalRawId', 'xtalEn', 'vSum8'])

myTree_noCh=tree_noCh.arrays(['runNumber', 'lumiBlock', 'eventNumber', 'xtalRawId', 'xtalEn', 'vSum8'])


run0=tree0.array('runNumber')
lumi0=tree0.array('lumiBlock')
evt0=tree0.array('eventNumber')
rawId0=zero_pad(tree0.array('xtalRawId'),5)
xtalEn0=zero_pad(tree0.array('xtalEn'),5)
vSum80=zero_pad(tree0.array('vSum8'), 5)

rawId0_0=rawId0[:,0]
rawId0_1=rawId0[:,1]
rawId0_2=rawId0[:,2]
rawId0_3=rawId0[:,3]
rawId0_4=rawId0[:,4]

xtalEn0_0=xtalEn0[:,0]
xtalEn0_1=xtalEn0[:,1]
xtalEn0_2=xtalEn0[:,2]
xtalEn0_3=xtalEn0[:,3]
xtalEn0_4=xtalEn0[:,4]

vSum80_0=vSum80[:,0]
vSum80_1=vSum80[:,1]
vSum80_2=vSum80[:,2]
vSum80_3=vSum80[:,3]
vSum80_4=vSum80[:,4]



run_noCh=tree_noCh.array('runNumber')
lumi_noCh=tree_noCh.array('lumiBlock')
evt_noCh=tree_noCh.array('eventNumber')
rawId_noCh=zero_pad(tree_noCh.array('xtalRawId'),5)
xtalEn_noCh=zero_pad(tree_noCh.array('xtalEn'),5)
vSum8_noCh=zero_pad(tree_noCh.array('vSum8'), 5)

rawId_noCh_0=rawId_noCh[:,0]
rawId_noCh_1=rawId_noCh[:,1]
rawId_noCh_2=rawId_noCh[:,2]
rawId_noCh_3=rawId_noCh[:,3]
rawId_noCh_4=rawId_noCh[:,4]

xtalEn_noCh_0=xtalEn_noCh[:,0]
xtalEn_noCh_1=xtalEn_noCh[:,1]
xtalEn_noCh_2=xtalEn_noCh[:,2]
xtalEn_noCh_3=xtalEn_noCh[:,3]
xtalEn_noCh_4=xtalEn_noCh[:,4]

vSum8_noCh_0=vSum8_noCh[:,0]
vSum8_noCh_1=vSum8_noCh[:,1]
vSum8_noCh_2=vSum8_noCh[:,2]
vSum8_noCh_3=vSum8_noCh[:,3]
vSum8_noCh_4=vSum8_noCh[:,4]



print len(xtalEn0), xtalEn0_1
#mylen=[]
#for ii in xrange(len(xtalEn0)):
#    mylen.append( len(xtalEn0[ii]))
#print max(mylen)
#print type(run0)

print run0, lumi0, evt0

df=pd.DataFrame(np.array([run0, lumi0, evt0, rawId0_0,rawId0_1, rawId0_2, rawId0_3, rawId0_4, xtalEn0_0, xtalEn0_1, xtalEn0_2, xtalEn0_3, xtalEn0_4, vSum80_0,  vSum80_1,  vSum80_2,  vSum80_3,  vSum80_4]).T, columns=['runNumber', 'lumiBlock', 'eventNumber', 'xtalRawId0', 'xtalRawId1', 'xtalRawId2', 'xtalRawId3', 'xtalRawId4','xtalEn0','xtalEn1','xtalEn2','xtalEn3','xtalEn4', 'vSum8_0', 'vSum8_1', 'vSum8_2', 'vSum8_3','vSum8_4' ] )

df_noCh=pd.DataFrame(np.array([run_noCh, lumi_noCh, evt_noCh, rawId_noCh_0,rawId_noCh_1, rawId_noCh_2, rawId_noCh_3, rawId_noCh_4, xtalEn_noCh_0, xtalEn_noCh_1, xtalEn_noCh_2, xtalEn_noCh_3, xtalEn_noCh_4,vSum8_noCh_0,vSum8_noCh_1,vSum8_noCh_2,vSum8_noCh_3,vSum8_noCh_4]).T, columns=['runNumber', 'lumiBlock', 'eventNumber', 'xtalRawId0', 'xtalRawId1', 'xtalRawId2', 'xtalRawId3', 'xtalRawId4','xtalEn0','xtalEn1','xtalEn2','xtalEn3','xtalEn4', 'vSum8_0', 'vSum8_1', 'vSum8_2', 'vSum8_3','vSum8_4'])

dfTot=pd.merge(df, df_noCh, how='inner', on=['runNumber', 'lumiBlock', 'eventNumber'], indicator=True)

dfTot= dfTot.loc[:, (dfTot != 0).any(axis=0)]
dfTot.rename(columns={'xtalRawId0_x':'xtalRawId0_recov', 'xtalRawId1_x':'xtalRawId1_recov', 'xtalRawId2_x':'xtalRawId2_recov', 'xtalRawId3_x':'xtalRawId3_recov', 'xtalRawId4_x':'xtalRawId4_recov', 'xtalEn0_x':'xtalEn0_recov', 'xtalEn1_x':'xtalEn1_recov', 'xtalEn2_x':'xtalEn2_recov', 'xtalEn3_x':'xtalEn3_recov', 'xtalEn4_x':'xtalEn4_recov', 'xtalRawId0_y':'xtalRawId0_true', 'xtalRawId1_y':'xtalRawId1_true', 'xtalRawId2_y':'xtalRawId2_true', 'xtalRawId3_y':'xtalRawId3_true', 'xtalRawId4_y':'xtalRawId4_true', 'xtalEn0_y':'xtalEn0_true', 'xtalEn1_y':'xtalEn1_true', 'xtalEn2_y':'xtalEn2_true', 'xtalEn3_y':'xtalEn3_true', 'xtalEn4_y':'xtalEn4_true', 'vSum8_0_x':'vSum8_0_recov','vSum8_1_x':'vSum8_1_recov','vSum8_2_x':'vSum8_2_recov','vSum8_3_x':'vSum8_3_recov','vSum8_4_x':'vSum8_4_recov', 'vSum8_0_y':'vSum8_0_true','vSum8_1_y':'vSum8_1_true','vSum8_2_y':'vSum8_2_true','vSum8_3_y':'vSum8_3_true','vSum8_4_y':'vSum8_4_true'  }, inplace=True)

#print dfTot.head()
#print dfTot.loc[:, ['xtalEn0_recov','xtalEn0_true','xtalEn1_recov','xtalEn1_true','xtalEn1_recov','xtalEn1_true', 'vSum8_0_true', 'vSum8_0_recov']]

dfDiff=dfTot[abs(dfTot.xtalEn0_recov-dfTot.xtalEn0_true)/dfTot.xtalEn0_true>0.50]
print dfDiff.shape[0]
print dfDiff.loc[:,  ['xtalEn0_recov','xtalEn0_true','xtalEn1_recov','xtalEn1_true','xtalEn1_recov','xtalEn1_true', 'vSum8_0_true', 'vSum8_0_recov']]

outfile=ROOT.TFile.Open("outfile_recovVsTrue.root", "RECREATE")
outfile.cd()

recov_vs_trueEn=ROOT.TH2F("recov_vs_trueEn", "recov_vs_trueEn", 100, 0, 500, 100, 0, 500)
diff_vs_sum8=ROOT.TH2F("diff_vs_sum8", "recov - true En vs sum8", 200, 0, 1000, 100, -250, 250)

true_vs_sum8=ROOT.TH2F("true_vs_sum8", "true_vs_sum8", 200, 0, 1000, 100, 0, 500)

for row in dfTot.itertuples():
    #print row.xtalEn0_true, row.xtalEn0_recov, row.xtalEn1_true, row.xtalEn1_recov,row.xtalEn2_true, row.xtalEn2_recov,row.xtalEn3_true,  row.xtalEn3_recov,row.xtalEn4_true,  row.xtalEn4_recov
    if row.xtalEn0_true!=0 or row.xtalEn0_recov!=0:
        recov_vs_trueEn.Fill(row.xtalEn0_true, row.xtalEn0_recov)
        diff_vs_sum8.Fill(row.vSum8_0_true, row.xtalEn0_recov-row.xtalEn0_true)
        true_vs_sum8.Fill(row.vSum8_0_true, row.xtalEn0_true)
    if row.xtalEn1_true!=0 and row.xtalEn1_recov!=0:
        recov_vs_trueEn.Fill(row.xtalEn1_true, row.xtalEn1_recov)
        diff_vs_sum8.Fill(row.vSum8_1_true, row.xtalEn1_recov-row.xtalEn1_true)
        true_vs_sum8.Fill(row.vSum8_1_true, row.xtalEn1_true)
    if row.xtalEn2_true!=0 and row.xtalEn2_recov!=0:
        recov_vs_trueEn.Fill(row.xtalEn2_true, row.xtalEn2_recov)
        diff_vs_sum8.Fill(row.vSum8_2_true, row.xtalEn2_recov-row.xtalEn2_true)
        true_vs_sum8.Fill(row.vSum8_2_true, row.xtalEn2_true)

    if row.xtalEn3_true!=0 and row.xtalEn3_recov!=0:
        recov_vs_trueEn.Fill(row.xtalEn3_true, row.xtalEn3_recov)
        diff_vs_sum8.Fill(row.vSum8_3_true, row.xtalEn3_recov-row.xtalEn3_true)
        true_vs_sum8.Fill(row.vSum8_3_true, row.xtalEn3_true)

    if row.xtalEn4_true!=0 and row.xtalEn4_recov!=0:
        recov_vs_trueEn.Fill(row.xtalEn4_true, row.xtalEn4_recov)
        diff_vs_sum8.Fill(row.vSum8_4_true, row.xtalEn4_recov-row.xtalEn4_true)
        true_vs_sum8.Fill(row.vSum8_4_true, row.xtalEn4_true)


outfile.cd()
recov_vs_trueEn.Write()
diff_vs_sum8.Write()
ROOT.gStyle.SetOptStat(0)

canvas = ROOT.TCanvas("c", "c", 600, 600)

canvas.Draw()
canvas.SetLeftMargin(0.15)
canvas.SetRightMargin(0.15)

recov_vs_trueEn.Draw("COLZ")
recov_vs_trueEn.GetXaxis().SetTitle("true En [GeV]")
recov_vs_trueEn.GetYaxis().SetTitle("recov En [GeV]")
canvas.SaveAs("plots/recov_vs_trueEn.png")
canvas.SaveAs("plots/recov_vs_trueEn.pdf")
c2=canvas.Clone()
c2.SetName("c_recov_vs_trueEn")
c2.Write()


canvas.Clear()

diff_vs_sum8.Draw("COLZ")
diff_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
diff_vs_sum8.GetYaxis().SetTitle("recov - true En [GeV]")
canvas.SaveAs("plots/diff_vs_sum8.pdf")
canvas.SaveAs("plots/diff_vs_sum8.png")
c3=canvas.Clone()
c3.SetName("c_diff_vs_sum8")
c3.Write()

canvas.Clear()

pfx = diff_vs_sum8.Rebin2D(2,2).ProfileX()
pfx.Draw("E")
pfx.GetXaxis().SetTitle("sum8 [GeV]")
pfx.GetYaxis().SetTitle("recov - true En [GeV]")
pfx.Write()
canvas.SaveAs("plots/profile_diff_vs_sum8.pdf")
canvas.SaveAs("plots/profile_diff_vs_sum8.png")
c4=canvas.Clone()
c4.SetName("c_profile_diff_vs_sum8")
c4.Write()


canvas.Clear()

true_vs_sum8.Draw("COLZ")
true_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
true_vs_sum8.GetYaxis().SetTitle("true En [GeV]")
canvas.SaveAs("plots/true_vs_sum8.png")
canvas.SaveAs("plots/true_vs_sum8.pdf")
c5=canvas.Clone()
c5.SetName("c_true_vs_sum8")
c5.Write()



outfile.Close()
