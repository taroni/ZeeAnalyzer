import ROOT
import uproot
import numpy as np
import pandas as pd
import array as arr

import argparse

file0=ROOT.TFile.Open("electronTreeZEE_PierreTag_sum8gt0_xtalInfo.root","READ")
#file0=ROOT.TFile.Open("electronTreeZEE_PierreTag_sum8gt20_xtalInfo.root", "READ")
#file0=ROOT.TFile.Open("electronTreeZEE_NoChange_sum8gt20_xtalInfo.root","READ")
tree0 =file0.Get("ntupler/selected")

outfile=ROOT.TFile("output_Neighbours_0GeV.root", "RECREATE")
#outfile=ROOT.TFile("output_Neighbours_20GeV.root", "RECREATE")
#outfile=ROOT.TFile("output_Neighbours_noChange.root", "RECREATE")

xtalEn=ROOT.TH1F("xtalEn", "xtalEn", 200, 0, 1000)
xtalEn_vs_sum8=ROOT.TH2F("xtalEn_vs_sum8", "xtalEn_vs_sum8", 100, 0, 1000, 100, 0, 500)
xtalEnOverTotEn=ROOT.TH1F("xtalEnOverTotEn", "xtalEnOverTotEn", 100, 0, 1)
xtalEn_vs_firstNeigh=ROOT.TH2F("xtalEn_vs_firstNeigh","xtalEn_vs_firstNeigh", 100, 0, 500, 100, 0, 500)
xtalEn_gtfirst = ROOT.TH1F("xtalEn_gtfirst","xtalEn_gtfirst", 100, 0, 500)
xtalEn_gtfirst_notCluster= ROOT.TH1F("xtalEn_gtfirst_notCluster","xtalEn_gtfirst_notCluster", 100, 0, 500)
nClust_xtalEngtfist_vs_nXtals=ROOT.TH2F("nClust_xtalEngtfist_vs_nXtals", "nClust_xtalEngtfist_vs_nXtals", 10, 0, 10, 10, 0, 10)
nClust_vs_nXtals = ROOT.TH2F("nClust_vs_nXtals", "nClust_vs_nXtals", 10, 0, 10, 10, 0, 10)

hSum8=ROOT.TH1F("hSum8", "sum8", 200, 0, 1000)
hAve8=ROOT.TH1F("hAve8", "ave8", 200, 0, 1000)
hE1overSum8=ROOT.TH1F("hE1overSum8", "E1overSum8", 100, 0., 1.)
hE2overSum8=ROOT.TH1F("hE2overSum8", "E2overSum8", 100, 0., 1.)
hE3overSum8=ROOT.TH1F("hE3overSum8", "E3overSum8", 100, 0., 1.)
hE4overSum8=ROOT.TH1F("hE4overSum8", "E4overSum8", 100, 0., 1.)
hE6overSum8=ROOT.TH1F("hE6overSum8", "E6overSum8", 100, 0., 1.)
hE7overSum8=ROOT.TH1F("hE7overSum8", "E7overSum8", 100, 0., 1.)
hE8overSum8=ROOT.TH1F("hE8overSum8", "E8overSum8", 100, 0., 1.)
hE9overSum8=ROOT.TH1F("hE9overSum8", "E9overSum8", 100, 0., 1.)

minTwoMaxTwoRatio=ROOT.TH1F("minTwoMaxTwoRatio", "minTwoMaxTwoRatio", 100, 0.,1.)
minTwoMaxTwoRatio_vs_sum8=ROOT.TH2F("minTwoMaxTwoRatio_vs_sum8", "minTwoMaxTwoRatio_vs_sum8", 100, 0, 1000, 100, 0, 1)
minTwoMaxTwoRatio_vs_invSum8=ROOT.TH2F("minTwoMaxTwoRatio_vs_invSum8", "minTwoMaxTwoRatio_vs_invSum8", 100, 0, 0.2, 100, 0, 1)
maxTwoMinTwoRatio_vs_sum8=ROOT.TH2F("maxTwoMinTwoRatio_vs_sum8", "maxTwoMinTwoRatio_vs_sum8", 100, 0, 1000, 100, 0, 100)
maxTwoMinTwoRatio_vs_xtalEn=ROOT.TH2F("maxTwoMinTwoRatio_vs_xtalEn", "maxTwoMinTwoRatio_vs_xtalEn", 100, 0, 1000, 100, 0, 100)

cornerOverCross=ROOT.TH1F("cornerOverCross", "cornerOverCross", 600, 0, 12)
cornerOverCross_vs_sum8=ROOT.TH2F("cornerOverCross_vs_sum8", "cornerOverCross_vs_sum8", 100, 0, 1000, 100, 0, 1)
cornerOverCross_vs_invSum8=ROOT.TH2F("cornerOverCross_vs_invSum8", "cornerOverCross_vs_invSum8", 100, 0, 0.2, 100, 0, 1)

firstNeigh=ROOT.TH1F("firstNeigh", "firstNeigh", 200, 0., 1000.)
secondNeigh=ROOT.TH1F("secondNeigh", "secondNeigh", 200, 0., 1000.)
thirdNeigh=ROOT.TH1F("thirdNeigh", "thirdNeigh", 200, 0., 1000.)
fourthNeigh=ROOT.TH1F("fourthNeigh", "fourthNeigh", 200, 0., 1000.)

firstEnFrac=ROOT.TH1F("firstEnFrac", "firstEnFrac", 100, 0., 1.)
seconEnFrac=ROOT.TH1F("secondEnFrac", "secondEnFrac", 100, 0., 1.)
thirdEnFrac=ROOT.TH1F("thirdEnFrac", "thirdEnFrac", 100, 0., 1.)
fourthEnFrac=ROOT.TH1F("fourthEnFrac", "fourthEnFrac", 100, 0., 1.)
fifthEnFrac=ROOT.TH1F("fifthEnFrac", "fifthEnFrac", 100, 0., 1.)
sixthEnFrac=ROOT.TH1F("sixthEnFrac", "sixthEnFrac", 100, 0., 1.)
seventhEnFrac=ROOT.TH1F("seventhEnFrac", "seventhEnFrac", 100, 0., 1.)
eighthEnFrac=ROOT.TH1F("eighthEnFrac", "eighthEnFrac", 100, 0., 1.)

first_overlast6=ROOT.TH1F("first_overlast6","first_overlast6", 200, 0, 10.)
firstsecond_overlast6=ROOT.TH1F("firstsecond_overlast6", "firstsecond_overlast6", 200, 0, 10.)
firstsecondEnFrac=ROOT.TH1F("firstsecondEnFrac", "firstsecondEnFrac", 100, 0, 1)
secondEnFrac_vs_firstEnFrac=ROOT.TH2F("secondEnFrac_vs_firstEnFrac", "secondEnFrac_vs_firstEnFrac", 100, 0, 1, 100, 0, 1)
thirdEnFrac_vs_firstEnFrac=ROOT.TH2F("thirdEnFrac_vs_firstEnFrac", "thirdEnFrac_vs_firstEnFrac", 100, 0, 1, 100, 0, 1)

                                  
secondfirstEnFrAsymmetry=ROOT.TH1F("secondfirstEnFrAsymmetry","secondfirstEnFrAsymmetry", 200, 0, 1)
thirdfirstEnFrAsymmetry=ROOT.TH1F("thirdfirstEnFrAsymmetry","thirdfirstEnFrAsymmetry", 200, 0, 1)
fourthfirstEnFrAsymmetry=ROOT.TH1F("fourthfirstEnFrAsymmetry","fourthfirstEnFrAsymmetry", 200, 0, 1)
fifthfirstEnFrAsymmetry=ROOT.TH1F("fifthfirstEnFrAsymmetry","fifthfirstEnFrAsymmetry", 200, 0, 1)
sixthfirstEnFrAsymmetry=ROOT.TH1F("sixthfirstEnFrAsymmetry","sixthfirstEnFrAsymmetry", 200, 0, 1)
seventhfirstEnFrAsymmetry=ROOT.TH1F("seventhfirstEnFrAsymmetry","seventhfirstEnFrAsymmetry", 200, 0, 1)
eighthfirstEnFrAsymmetry=ROOT.TH1F("eighthfirstEnFrAsymmetry","eighthfirstEnFrAsymmetry", 200, 0, 1)

secondfirstEnAsymmetry_vs_sum8 =ROOT.TH2F("secondfirstEnAsymmetry_vs_sum8","secondfirstEnAsymmetry_vs_sum8"  , 200, 0, 1000, 200, 0, 1)
thirdfirstEnAsymmetry_vs_sum8  =ROOT.TH2F("thirdfirstEnAsymmetry_vs_sum8","thirdfirstEnAsymmetry_vs_sum8"    , 200, 0, 1000, 200, 0, 1)
fourthfirstEnAsymmetry_vs_sum8 =ROOT.TH2F("fourthfirstEnAsymmetry_vs_sum8","fourthfirstEnAsymmetry_vs_sum8"  , 200, 0, 1000, 200, 0, 1)
fifthfirstEnAsymmetry_vs_sum8  =ROOT.TH2F("fifthfirstEnAsymmetry_vs_sum8","fifthfirstEnAsymmetry_vs_sum8"    , 200, 0, 1000, 200, 0, 1)
sixthfirstEnAsymmetry_vs_sum8  =ROOT.TH2F("sixthfirstEnAsymmetry_vs_sum8","sixthfirstEnAsymmetry_vs_sum8"    , 200, 0, 1000, 200, 0, 1)
seventhfirstEnAsymmetry_vs_sum8=ROOT.TH2F("seventhfirstEnAsymmetry_vs_sum8","seventhfirstEnAsymmetry_vs_sum8", 200, 0, 1000, 200, 0, 1)
eighthfirstEnAsymmetry_vs_sum8 =ROOT.TH2F("eighthfirstEnAsymmetry_vs_sum8","eighthfirstEnAsymmetry_vs_sum8"  , 200, 0, 1000, 200, 0, 1)


second_vs_firstNeighIndex=ROOT.TH2F("second_vs_firstNeighIndex", "second_vs_firstNeighIndex", 10, 0, 10, 10, 0, 10)


firstNeighIndex=ROOT.TH1F("firstNeighIndex", "firstNeighIndex", 10, 0., 10.)
secondNeighIndex=ROOT.TH1F("secondNeighIndex", "secondNeighIndex", 10, 0., 10.)

firstNeighIndex_vs_xtalEn=ROOT.TH2F("firstNeighIndex_vs_xtalEn", "firstNeighIndex_vs_xtalEn", 100, 0, 1000, 10, 0., 10.)
secondNeighIndex_vs_xtalEn=ROOT.TH2F("secondNeighIndex_vs_xtalEn", "secondNeighIndex_vs_xtalEn", 100, 0, 1000, 10, 0., 10.)

firstNeighIndex_vs_cornerOverCross=ROOT.TH2F("firstNeighIndex_vs_cornerOverCross", "firstNeighIndex_vs_cornerOverCross",150, 0, 1.5, 10, 0, 10)
secondNeighIndex_vs_cornerOverCross=ROOT.TH2F("secondNeighIndex_vs_cornerOverCross", "secondNeighIndex_vs_cornerOverCross", 150, 0, 1.5,  10, 0, 10)



second_vs_first=ROOT.TH2F("second_vs_first", "second_vs_first", 100, 0, 1000, 50, 0, 500)
third_vs_first=ROOT.TH2F("third_vs_first", "third_vs_first", 100, 0, 1000, 50, 0, 500)

first_vs_sum8=ROOT.TH2F("first_vs_sum8", "first_vs_sum8",  100, 0,1000,100, 0,1000)
second_vs_sum8=ROOT.TH2F("second_vs_sum8", "second_vs_sum8",  100, 0,1000,100, 0,1000)
third_vs_sum8=ROOT.TH2F("third_vs_sum8", "third_vs_sum8",  100, 0,1000,100, 0,1000)

hE1oS_vs_sum8=ROOT.TH2F("hE1oS_vs_sum8", "hE1oS_vs_sum8", 100, 0,1000, 100, 0,1)
hE2oS_vs_sum8=ROOT.TH2F("hE2oS_vs_sum8", "hE2oS_vs_sum8", 100, 0,1000, 100, 0,1)
hE3oS_vs_sum8=ROOT.TH2F("hE3oS_vs_sum8", "hE3oS_vs_sum8", 100, 0,1000, 100, 0,1)
hE4oS_vs_sum8=ROOT.TH2F("hE4oS_vs_sum8", "hE4oS_vs_sum8", 100, 0,1000, 100, 0,1)
hE6oS_vs_sum8=ROOT.TH2F("hE6oS_vs_sum8", "hE6oS_vs_sum8", 100, 0,1000, 100, 0,1)
hE7oS_vs_sum8=ROOT.TH2F("hE7oS_vs_sum8", "hE7oS_vs_sum8", 100, 0,1000, 100, 0,1)
hE8oS_vs_sum8=ROOT.TH2F("hE8oS_vs_sum8", "hE8oS_vs_sum8", 100, 0,1000, 100, 0,1)
hE9oS_vs_sum8=ROOT.TH2F("hE9oS_vs_sum8", "hE9oS_vs_sum8", 100, 0,1000, 100, 0,1)


nXtal90=ROOT.TH1F("nXtal90", "nXtal90", 9, 0, 9); 
nXtal80=ROOT.TH1F("nXtal80", "nXtal80", 9, 0, 9); 
nXtal70=ROOT.TH1F("nXtal70", "nXtal70", 9, 0, 9); 
nXtal50=ROOT.TH1F("nXtal50", "nXtal50", 9, 0, 9); 
nXtal30=ROOT.TH1F("nXtal30", "nXtal30", 9, 0, 9); 

leftOverSum8=ROOT.TH1F("leftOverSum8", "leftOverSum8", 100, 0, 1)
rightOverSum8=ROOT.TH1F("rightOverSum8", "rightOverSum8", 100, 0, 1)
leftRightOverSum8=ROOT.TH1F("leftRightOverSum8", "leftRightOverSum8", 100, 0, 1)
leftRightOverSum8_vs_sum8=ROOT.TH2F("leftRightOverSum8_vs_sum8", "leftRightOverSum8_vs_sum8", 100, 0, 1000, 100, 0, 1)
left_vs_sum8=ROOT.TH2F("left_vs_sum8", "left_vs_sum8", 100, 0, 1000, 50, 0,500) 
right_vs_sum8=ROOT.TH2F("right_vs_sum8", "right_vs_sum8", 100, 0, 1000, 50, 0,500) 
leftRight_vs_sum8=ROOT.TH2F("leftRight_vs_sum8", "leftRight_vs_sum8", 100, 0, 1000, 50, 0,500) 

upOverSum8=ROOT.TH1F("upOverSum8", "upOverSum8", 100, 0, 1)
downOverSum8=ROOT.TH1F("downOverSum8", "downOverSum8", 100, 0, 1)
upDownOverSum8=ROOT.TH1F("upDownOverSum8", "upDownOverSum8", 100, 0, 1)
upDownOverSum8_vs_sum8=ROOT.TH2F("upDownOverSum8_vs_sum8", "upDownOverSum8_vs_sum8", 100, 0, 1000, 100, 0, 1)
up_vs_sum8=ROOT.TH2F("up_vs_sum8", "up_vs_sum8", 100, 0, 1000, 50, 0,500) 
down_vs_sum8=ROOT.TH2F("down_vs_sum8", "down_vs_sum8", 100, 0, 1000, 50, 0,500) 
upDown_vs_sum8=ROOT.TH2F("upDown_vs_sum8", "upDown_vs_sum8", 100, 0, 1000, 50, 0,500) 



eleEn=ROOT.TH1F("eleEn", "eleEn", 150, 0, 1500)
eleRawEn=ROOT.TH1F("eleRawEn", "eleRawEn", 150, 0, 1500)
eleEn_leftRightOverSum8LT0p3=ROOT.TH1F("eleEn_leftRightOverSum8LT0p3", "eleEn_leftRightOverSum8LT0p3", 150, 0, 1500)
eleRawEn_leftRightOverSum8LT0p3=ROOT.TH1F("eleRawEn_leftRightOverSum8LT0p3", "eleRawEn_leftRightOverSum8LT0p3", 150, 0, 1500)
eleEn_leftRightOverSum8LT0p5=ROOT.TH1F("eleEn_leftRightOverSum8LT0p5", "eleEn_leftRightOverSum8LT0p5", 150, 0, 1500)
eleRawEn_leftRightOverSum8LT0p5=ROOT.TH1F("eleRawEn_leftRightOverSum8LT0p5", "eleRawEn_leftRightOverSum8LT0p5", 150, 0, 1500)

for row in tree0:
    
    e1Recov=row.e1IsRecovered
    e2Recov=row.e2IsRecovered

    e1En=row.e1Energy
    e2En=row.e1Energy

    e1RawEn=row.e1RawEnergy
    e2RawEn=row.e2RawEnergy
    #print 'number of recovered xtals %i, e1Recov %r, e2Recov %r' %(len(row.xtalEn), e1Recov, e2Recov)

    for ixtal in range(len(row.xtalEn)):
        xtalEn.Fill(row.xtalEn[ixtal])
        xtalEn_vs_sum8.Fill(row.vSum8[ixtal],row.xtalEn[ixtal])
        
    
        energies=[]
        energies.extend([row.vNeigh0[ixtal],row.vNeigh1[ixtal],row.vNeigh2[ixtal],row.vNeigh3[ixtal],row.vNeigh5[ixtal],row.vNeigh6[ixtal],row.vNeigh7[ixtal],row.vNeigh8[ixtal]])
        
        leftEn=row.vNeigh0[ixtal]+row.vNeigh1[ixtal]+row.vNeigh2[ixtal]
        rightEn=row.vNeigh6[ixtal]+row.vNeigh7[ixtal]+row.vNeigh8[ixtal]

        downEn=row.vNeigh0[ixtal]+row.vNeigh3[ixtal]+row.vNeigh6[ixtal]
        upEn=row.vNeigh2[ixtal]+row.vNeigh5[ixtal]+row.vNeigh8[ixtal]

        maxIndex=energies.index(max(energies))+1
        if maxIndex>=5: maxIndex+=1
        firstNeighIndex.Fill(maxIndex)
        copiedEnergies=[x for x in energies if x!=max(energies)]
        #print energies
        #print copiedEnergies
        secondmax=max(copiedEnergies)
        #print max(energies), secondmax
        secondMaxIndex=energies.index(secondmax)+1
        if secondMaxIndex>=5: secondMaxIndex+=1
        secondNeighIndex.Fill(secondMaxIndex)
        #print maxIndex, secondMaxIndex
        second_vs_firstNeighIndex.Fill(secondMaxIndex, maxIndex)

        orderedEn=[x for x in energies]
        orderedEn.sort(reverse=True)
        #if len(row.xtalEn)>1 : print 'number of recovered xtals %i, energy %f, %f, %f, %f, %f' %(len(row.xtalEn),row.xtalEn[ixtal],orderedEn[0], orderedEn[1], orderedEn[2], orderedEn[3])
        #print orderedEn
        firstNeigh.Fill(orderedEn[0])
        secondNeigh.Fill(orderedEn[1])
        second_vs_first.Fill(orderedEn[0], orderedEn[1])
        
        thirdNeigh.Fill(orderedEn[2]) 
        third_vs_first.Fill(orderedEn[0], orderedEn[2])
        fourthNeigh.Fill(orderedEn[3]) 
               
        firstNeighIndex_vs_xtalEn.Fill(row.xtalEn[ixtal],maxIndex)
        secondNeighIndex_vs_xtalEn.Fill(row.xtalEn[ixtal],secondMaxIndex)

        xtalEnOverTotEn.Fill(row.xtalEn[ixtal]/(row.xtalEn[ixtal]+row.vSum8[ixtal]))
        xtalEn_vs_firstNeigh.Fill(orderedEn[0], row.xtalEn[ixtal])
        tmpXtalEn=list(set([row.xtalEn[ii] for ii in xrange(0, len(row.xtalEn))]))
        
        nClus=0
        tmpEn=0
        xtalList=[x for x in row.xtalEn]
        ii=0
        while ii < len(xtalList):
            #print  x, tmpEn, xtalList, nClus
            if tmpEn==x:
                xtalList = [ elem for elem in xtalList if elem!=x]
                nClus+=1
            tmpEn=x
            ii+=1
            
        nClust_vs_nXtals.Fill(len(row.xtalEn), nClus )

        if len(row.xtalEn)>1: print ixtal, len(row.xtalEn), len(tmpXtalEn), row.xtalEn[ixtal],  tmpXtalEn, orderedEn[0], row.xtalIeta[ixtal], row.xtalIphi[ixtal]
        if row.xtalEn[ixtal]>orderedEn[0]: 
            xtalEn_gtfirst.Fill(row.xtalEn[ixtal])
            if row.xtalEn[ixtal] in xtalList:
                xtalEn_gtfirst_notCluster.Fill(row.xtalEn[ixtal])

            #if len(row.xtalEn)!=1 and len(xtalList)<len(row.xtalEn) :
            #    nClust_xtalEngtfist_vs_nXtals.Fill(len(row.xtalEn), nClus)
            


        first_vs_sum8.Fill(row.vSum8[ixtal], orderedEn[0])
        second_vs_sum8.Fill(row.vSum8[ixtal],orderedEn[1])
        third_vs_sum8.Fill(row.vSum8[ixtal],orderedEn[2])

        firstEnFrac.Fill( orderedEn[0]/row.vSum8[ixtal])
        seconEnFrac.Fill( orderedEn[1]/row.vSum8[ixtal])
        thirdEnFrac.Fill( orderedEn[2]/row.vSum8[ixtal])
        fourthEnFrac.Fill( orderedEn[3]/row.vSum8[ixtal])
        fifthEnFrac.Fill( orderedEn[4]/row.vSum8[ixtal])
        sixthEnFrac.Fill( orderedEn[5]/row.vSum8[ixtal])
        seventhEnFrac.Fill( orderedEn[6]/row.vSum8[ixtal])
        eighthEnFrac.Fill( orderedEn[7]/row.vSum8[ixtal])

        first_overlast6.Fill(orderedEn[0]/(orderedEn[2]+orderedEn[2]+orderedEn[4]+orderedEn[5]+orderedEn[6]+orderedEn[7]))
        firstsecond_overlast6.Fill((orderedEn[0]+orderedEn[1])/(orderedEn[2]+orderedEn[2]+orderedEn[4]+orderedEn[5]+orderedEn[6]+orderedEn[7]))
        firstsecondEnFrac.Fill((orderedEn[0]+orderedEn[1])/row.vSum8[ixtal])
        secondEnFrac_vs_firstEnFrac.Fill(orderedEn[0]/row.vSum8[ixtal], orderedEn[1]/row.vSum8[ixtal])
        thirdEnFrac_vs_firstEnFrac.Fill(orderedEn[0]/row.vSum8[ixtal], orderedEn[2]/row.vSum8[ixtal])
        secondfirstEnFrAsymmetry.Fill((orderedEn[0]-orderedEn[1])/(orderedEn[0]+orderedEn[1]))
        thirdfirstEnFrAsymmetry.Fill((orderedEn[0]-orderedEn[2])/(orderedEn[0]+orderedEn[2]))
        fourthfirstEnFrAsymmetry.Fill((orderedEn[0]-orderedEn[3])/(orderedEn[0]+orderedEn[3]))
        fifthfirstEnFrAsymmetry.Fill((orderedEn[0]-orderedEn[4])/(orderedEn[0]+orderedEn[4]))
        sixthfirstEnFrAsymmetry.Fill((orderedEn[0]-orderedEn[5])/(orderedEn[0]+orderedEn[5]))
        seventhfirstEnFrAsymmetry.Fill((orderedEn[0]-orderedEn[6])/(orderedEn[0]+orderedEn[6]))
        eighthfirstEnFrAsymmetry.Fill((orderedEn[0]-orderedEn[7])/(orderedEn[0]+orderedEn[7]))

        secondfirstEnAsymmetry_vs_sum8.Fill(row.vSum8[ixtal],(orderedEn[0]-orderedEn[1])/(orderedEn[0]+orderedEn[1]))
        thirdfirstEnAsymmetry_vs_sum8.Fill(row.vSum8[ixtal],(orderedEn[0]-orderedEn[2])/(orderedEn[0]+orderedEn[2]))
        fourthfirstEnAsymmetry_vs_sum8.Fill(row.vSum8[ixtal],(orderedEn[0]-orderedEn[3])/(orderedEn[0]+orderedEn[3]))
        fifthfirstEnAsymmetry_vs_sum8.Fill(row.vSum8[ixtal],(orderedEn[0]-orderedEn[4])/(orderedEn[0]+orderedEn[4]))
        sixthfirstEnAsymmetry_vs_sum8.Fill(row.vSum8[ixtal],(orderedEn[0]-orderedEn[5])/(orderedEn[0]+orderedEn[5]))
        seventhfirstEnAsymmetry_vs_sum8.Fill(row.vSum8[ixtal],(orderedEn[0]-orderedEn[6])/(orderedEn[0]+orderedEn[6]))
        eighthfirstEnAsymmetry_vs_sum8.Fill(row.vSum8[ixtal],(orderedEn[0]-orderedEn[7])/(orderedEn[0]+orderedEn[7]))


        maxTwo=orderedEn[0]+orderedEn[1]
        minTwo=orderedEn[-1]+orderedEn[-2]
        
        if maxTwo !=0:
            minTwoMaxTwoRatio.Fill(minTwo/maxTwo)
            minTwoMaxTwoRatio_vs_sum8.Fill(row.vSum8[ixtal],minTwo/maxTwo)
            minTwoMaxTwoRatio_vs_invSum8.Fill(1./row.vSum8[ixtal],minTwo/maxTwo)

        if minTwo!=0:
            maxTwoMinTwoRatio_vs_sum8.Fill(row.vSum8[ixtal], maxTwo/minTwo)
            maxTwoMinTwoRatio_vs_xtalEn.Fill(row.xtalEn[ixtal], maxTwo/minTwo)

        #print energies
        #print orderedEn
        corners=energies[0]+energies[2]+energies[5]+energies[7]
        cross=energies[1]+energies[3]+energies[4]+energies[6]
        #print corners, cross, row.vSum8[ixtal]
        if cross !=0:
            cornerOverCross.Fill(corners/cross)
            cornerOverCross_vs_sum8.Fill(row.vSum8[ixtal],corners/cross)

            cornerOverCross_vs_invSum8.Fill(1./row.vSum8[ixtal],corners/cross)
            

        firstNeighIndex_vs_cornerOverCross.Fill(corners/cross, maxIndex)
        secondNeighIndex_vs_cornerOverCross.Fill(corners/cross,secondMaxIndex)


        icount90=1
        icount80=1
        icount70=1
        icount50=1
        icount30=1
        mysum=0.
        for en in orderedEn:
            mysum+=en
            #print ixtal, mysum, row.vSum8[ixtal], icount50, icount70, icount80, icount90
            if mysum>=0.9*row.vSum8[ixtal]: continue
            icount90+=1
            if mysum>=0.8*row.vSum8[ixtal]: continue
            icount80+=1
            if mysum>=0.7*row.vSum8[ixtal]: continue
            icount70+=1
            if mysum>=0.5*row.vSum8[ixtal]: continue
            icount50+=1
            if mysum>=0.3*row.vSum8[ixtal]: continue
            icount30+=1

        #print icount30, icount50, icount70, icount80, icount90

        nXtal90.Fill(icount90)
        nXtal80.Fill(icount80)
        nXtal70.Fill(icount70)
        nXtal50.Fill(icount50)
        nXtal30.Fill(icount30)
        
        hSum8.Fill(row.vSum8[ixtal])
        hAve8.Fill(row.vAve[ixtal])
        hE1overSum8.Fill(row.vNeigh0[ixtal]/row.vSum8[ixtal])
        hE2overSum8.Fill(row.vNeigh1[ixtal]/row.vSum8[ixtal])
        hE3overSum8.Fill(row.vNeigh2[ixtal]/row.vSum8[ixtal])
        hE4overSum8.Fill(row.vNeigh3[ixtal]/row.vSum8[ixtal])
        hE6overSum8.Fill(row.vNeigh5[ixtal]/row.vSum8[ixtal])
        hE7overSum8.Fill(row.vNeigh6[ixtal]/row.vSum8[ixtal])
        hE8overSum8.Fill(row.vNeigh7[ixtal]/row.vSum8[ixtal])
        hE9overSum8.Fill(row.vNeigh8[ixtal]/row.vSum8[ixtal])


        hE1oS_vs_sum8.Fill(row.vSum8[ixtal],row.vNeigh0[ixtal]/row.vSum8[ixtal])
        hE2oS_vs_sum8.Fill(row.vSum8[ixtal],row.vNeigh1[ixtal]/row.vSum8[ixtal])
        hE3oS_vs_sum8.Fill(row.vSum8[ixtal],row.vNeigh2[ixtal]/row.vSum8[ixtal])
        hE4oS_vs_sum8.Fill(row.vSum8[ixtal],row.vNeigh3[ixtal]/row.vSum8[ixtal])
        hE6oS_vs_sum8.Fill(row.vSum8[ixtal],row.vNeigh5[ixtal]/row.vSum8[ixtal])
        hE7oS_vs_sum8.Fill(row.vSum8[ixtal],row.vNeigh6[ixtal]/row.vSum8[ixtal])
        hE8oS_vs_sum8.Fill(row.vSum8[ixtal],row.vNeigh7[ixtal]/row.vSum8[ixtal])
        hE9oS_vs_sum8.Fill(row.vSum8[ixtal],row.vNeigh8[ixtal]/row.vSum8[ixtal])

        leftOverSum8.Fill(leftEn/row.vSum8[ixtal])
        rightOverSum8.Fill(rightEn/row.vSum8[ixtal])
        leftRightOverSum8.Fill((leftEn+rightEn)/row.vSum8[ixtal])
        leftRightOverSum8_vs_sum8.Fill(row.vSum8[ixtal],(leftEn+rightEn)/row.vSum8[ixtal])
        
        left_vs_sum8.Fill(row.vSum8[ixtal], leftEn)
        right_vs_sum8.Fill(row.vSum8[ixtal],rightEn)
        leftRight_vs_sum8.Fill(row.vSum8[ixtal], leftEn+rightEn)

        upOverSum8.Fill(upEn/row.vSum8[ixtal])
        downOverSum8.Fill(downEn/row.vSum8[ixtal])
        upDownOverSum8.Fill((upEn+downEn)/row.vSum8[ixtal])
        upDownOverSum8_vs_sum8.Fill(row.vSum8[ixtal],(upEn+downEn)/row.vSum8[ixtal])
        up_vs_sum8.Fill(row.vSum8[ixtal], upEn)
        down_vs_sum8.Fill(row.vSum8[ixtal],downEn)
        upDown_vs_sum8.Fill(row.vSum8[ixtal], upEn+downEn)


        if e1Recov :
            eleEn.Fill(e1En)
            eleRawEn.Fill(e1RawEn)
            if (leftEn+rightEn)/row.vSum8[ixtal]<0.3:
                eleEn_leftRightOverSum8LT0p3.Fill(e1En)
                eleRawEn_leftRightOverSum8LT0p3.Fill(e1RawEn)
            if (leftEn+rightEn)/row.vSum8[ixtal]<0.5:
                eleEn_leftRightOverSum8LT0p5.Fill(e1En)
                eleRawEn_leftRightOverSum8LT0p5.Fill(e1RawEn)
        if e2Recov:
            eleEn.Fill(e2En)
            eleRawEn.Fill(e2RawEn)
            if (leftEn+rightEn)/row.vSum8[ixtal]<0.3:
                eleEn_leftRightOverSum8LT0p3.Fill(e2En)
                eleRawEn_leftRightOverSum8LT0p3.Fill(e2RawEn)
            if (leftEn+rightEn)/row.vSum8[ixtal]<0.5:
                eleEn_leftRightOverSum8LT0p5.Fill(e2En)
                eleRawEn_leftRightOverSum8LT0p5.Fill(e2RawEn)
  
            
            

outfile.cd()
xtalEn.Write()
hSum8.Write()
hAve8.Write()
hE1overSum8.Write()
hE2overSum8.Write()
hE3overSum8.Write()
hE4overSum8.Write()
hE6overSum8.Write()
hE7overSum8.Write()
hE8overSum8.Write()
hE9overSum8.Write()

firstNeigh.Write()
secondNeigh.Write()
thirdNeigh.Write()
fourthNeigh.Write()

firstNeighIndex.Write()
secondNeighIndex.Write()

second_vs_first.Write()
third_vs_first.Write()

first_vs_sum8.Write()
second_vs_sum8.Write()
third_vs_sum8.Write()

firstEnFrac.Write()
seconEnFrac.Write()
thirdEnFrac.Write()
fourthEnFrac.Write()
fifthEnFrac.Write()
sixthEnFrac.Write()
seventhEnFrac.Write()
eighthEnFrac.Write()

first_overlast6.Write()
firstsecond_overlast6.Write()
firstsecondEnFrac.Write()
secondEnFrac_vs_firstEnFrac.Write()
thirdEnFrac_vs_firstEnFrac.Write()
secondfirstEnFrAsymmetry.Write()
thirdfirstEnFrAsymmetry.Write()
fourthfirstEnFrAsymmetry.Write()
fifthfirstEnFrAsymmetry.Write()
sixthfirstEnFrAsymmetry.Write()
seventhfirstEnFrAsymmetry.Write()
eighthfirstEnFrAsymmetry.Write()
                                   
secondfirstEnAsymmetry_vs_sum8.Write()
thirdfirstEnAsymmetry_vs_sum8.Write()  
fourthfirstEnAsymmetry_vs_sum8.Write() 
fifthfirstEnAsymmetry_vs_sum8.Write()  
sixthfirstEnAsymmetry_vs_sum8.Write()  
seventhfirstEnAsymmetry_vs_sum8.Write()
eighthfirstEnAsymmetry_vs_sum8.Write() 

nXtal90.Write()
nXtal80.Write()
nXtal70.Write()
nXtal50.Write()
nXtal30.Write()

hE1oS_vs_sum8.Write()
hE2oS_vs_sum8.Write()
hE3oS_vs_sum8.Write()
hE4oS_vs_sum8.Write()
hE6oS_vs_sum8.Write()
hE7oS_vs_sum8.Write()
hE8oS_vs_sum8.Write()
hE9oS_vs_sum8.Write()

cornerOverCross.Write()
cornerOverCross_vs_sum8.Write()
minTwoMaxTwoRatio.Write()
minTwoMaxTwoRatio_vs_sum8.Write()
minTwoMaxTwoRatio_vs_invSum8.Write()
maxTwoMinTwoRatio_vs_sum8.Write()

firstNeighIndex_vs_xtalEn.Write()
secondNeighIndex_vs_xtalEn.Write()
maxTwoMinTwoRatio_vs_xtalEn.Write()

firstNeighIndex_vs_cornerOverCross.Write()
secondNeighIndex_vs_cornerOverCross.Write()

second_vs_firstNeighIndex.Write()


cornerOverCross_vs_invSum8.Write()

leftOverSum8.Write()
rightOverSum8.Write()
leftRightOverSum8.Write()
leftRightOverSum8_vs_sum8.Write()
left_vs_sum8.Write()
right_vs_sum8.Write()
leftRight_vs_sum8.Write()

upOverSum8.Write()
downOverSum8.Write()
upDownOverSum8.Write()
upDownOverSum8_vs_sum8.Write()
up_vs_sum8.Write()
down_vs_sum8.Write()
upDown_vs_sum8.Write()


eleEn.Write()
eleRawEn.Write()
eleEn_leftRightOverSum8LT0p3.Write()
eleRawEn_leftRightOverSum8LT0p3.Write()
eleEn_leftRightOverSum8LT0p5.Write()
eleRawEn_leftRightOverSum8LT0p5.Write()

xtalEn_vs_sum8.Write()

xtalEnOverTotEn.Write()
xtalEn_vs_firstNeigh.Write()
xtalEn_gtfirst.Write()
nClust_xtalEngtfist_vs_nXtals.Write()
nClust_vs_nXtals.Write()
xtalEn_gtfirst_notCluster.Write()

outfile.Close()
