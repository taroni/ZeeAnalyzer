import ROOT

infile = ROOT.TFile.Open("electronTreeZEE_PierreTag_official_xtalInfo.root","READ")

tree = infile.Get("ntupler/selected")

hRecovXtalSize = ROOT.TH1F('hRecovXtalsSize', 'hRecovXtalSize',100, 0, 100)
hGood_recov1                 = ROOT.TH1F("hGood_recov1", "hGood_recov1", 2, 0, 2)
hPoorReco_recov1             = ROOT.TH1F("hPoorReco_recov1", "hPoorReco_recov1", 2, 0, 2)
hOutOfTimE_recov1            = ROOT.TH1F("hOutOfTimE_recov1", "hOutOfTimE_recov1", 2, 0, 2)
hFaultyHardware_recov1       = ROOT.TH1F("hFaultyHardware_recov1", "hFaultyHardware_recov1", 2, 0, 2)
hNoisy_recov1                = ROOT.TH1F("hNoisy_recov1", "hNoisy_recov1", 2, 0, 2)
hPoorCalib_recov1            = ROOT.TH1F("hPoorCalib_recov1", "hPoorCalib_recov1", 2, 0, 2)
hSaturated_recov1            = ROOT.TH1F("hSaturated_recov1", "hSaturated_recov1", 2, 0, 2)
hLeadingEdgeRecovered_recov1 = ROOT.TH1F("hLeadingEdgeRecovered_recov1", "hLeadingEdgeRecovered_recov1", 2, 0, 2)
hNeighboursRecovered_recov1  = ROOT.TH1F("hNeighboursRecovered_recov1", "hNeighboursRecovered_recov1", 2, 0, 2)
hTowerRecovered_recov1       = ROOT.TH1F("hTowerRecovered_recov1", "hTowerRecovered_recov1", 2, 0, 2)
hDead_recov1                 = ROOT.TH1F("hDead_recov1", "hDead_recov1", 2, 0, 2)
hKilled_recov1               = ROOT.TH1F("hKilled_recov1", "hKilled_recov1", 2, 0, 2)
hTPSaturated_recov1          = ROOT.TH1F("hTPSaturated_recov1", "hTPSaturated_recov1", 2, 0, 2)
hL1SpikeFlag_recov1          = ROOT.TH1F("hL1SpikeFlag_recov1", "hL1SpikeFlag_recov1", 2, 0, 2)
hWeird_recov1                = ROOT.TH1F("hWeird_recov1", "hWeird_recov1", 2, 0, 2)
hDiWeird_recov1              = ROOT.TH1F("hDiWeird_recov1", "hDiWeird_recov1", 2, 0, 2)
hHasSwitchToGain6_recov1     = ROOT.TH1F("hHasSwitchToGain6_recov1", "hHasSwitchToGain6_recov1", 2, 0, 2)
hHasSwitchToGain1_recov1     = ROOT.TH1F("hHasSwitchToGain1_recov1", "hHasSwitchToGain1_recov1", 2, 0, 2)

#kGood
#kPoorReco         
#kOutOfTimE         
#kFaultyHardware      
#kNoisy               
#kPoorCalib          
#kSaturated           
#kLeadingEdgeRecovered
#kNeighboursRecovered 
#kTowerRecovered      
#kDead                
#kKilled              
#kTPSaturated         
#kL1SpikeFlag         
#kWeird               
#kDiWeird             
#kHasSwitchToGain6    
#kHasSwitchToGain1    
#kUnknown              




for evt in tree:
    
    kRecov=evt.kNeighboursRecovered
    hRecovXtalSize.Fill(len(kRecov))

    if len(kRecov) <=10:
        for ixtal in range(0, len(kRecov)):
            ##hGood_recov1.Fill(evt.kGood[ixtal])            
            ##hPoorReco_recov1.Fill(evt.kPoorReco[ixtal])            
            ##hOutOfTimE_recov1.Fill(evt.kOutOfTimE[ixtal])
            ##hFaultyHardware_recov1.Fill(evt.kFaultyHardware[ixtal])
            ##hNoisy_recov1.Fill(evt.kNoisy[ixtal])
            ##hPoorCalib_recov1.Fill(evt.kPoorCalib[ixtal])
            ##hSaturated_recov1.Fill(evt.kSaturated[ixtal])           
            ##hLeadingEdgeRecovered_recov1.Fill(evt.kLeadingEdgeRecovered[ixtal])
            ##hNeighboursRecovered_recov1.Fill(evt.kNeighboursRecovered[ixtal])
            ##hTowerRecovered_recov1.Fill(evt.kTowerRecovered[ixtal])      
            ##hDead_recov1.Fill(evt.kDead[ixtal])                
            ##hKilled_recov1.Fill(evt.kKilled[ixtal])              
            ##hTPSaturated_recov1.Fill(evt.kTPSaturated[ixtal])
            ##hL1SpikeFlag_recov1.Fill(evt.kL1SpikeFlag[ixtal])         
            ##hWeird_recov1.Fill(evt.kWeird[ixtal])
            ##hDiWeird_recov1.Fill(evt.kDiWeird[ixtal])             
            ##hHasSwitchToGain6_recov1.Fill(evt.kHasSwitchToGain6[ixtal])   
            ##hHasSwitchToGain1_recov1.Fill(evt.kHasSwitchToGain1[ixtal])     
    
            ##if evt.kGood[ixtal]==True or evt.kDead[ixtal]==True:
            ##print 'kGood %d, kPoorReco %d, kOutOfTimE %d, kFaultyHardware %d, kNoisy %d, kPoorCalib %d, kSaturated %d, kLeadingEdgeRecovered %d, kNeighboursRecovered %d, kTowerRecovered %d, kDead %d, kKilled %d, kTPSaturated %d,kL1SpikeFlag %d, kWeird %d, kDiWeird %d, kHasSwitchToGain6 %d, kHasSwitchToGain1 %d, kUnknown %d' %(   evt.kGood[ixtal],  evt.kPoorReco[ixtal],          evt.kOutOfTimE[ixtal],           evt.kFaultyHardware[ixtal],        evt.kNoisy[ixtal],                 evt.kPoorCalib[ixtal],            evt.kSaturated[ixtal],             evt.kLeadingEdgeRecovered[ixtal],  evt.kNeighboursRecovered[ixtal],   evt.kTowerRecovered[ixtal],        evt.kDead[ixtal],                 evt.kKilled[ixtal],                evt.kTPSaturated[ixtal],           evt.kL1SpikeFlag[ixtal],           evt.kWeird[ixtal],                 evt.kDiWeird[ixtal],               evt.kHasSwitchToGain6[ixtal],      evt.kHasSwitchToGain1[ixtal],      evt.kUnknown[ixtal] )
            print 'kGood %s, kPoorReco %s, kOutOfTimE %s, kFaultyHardware %s, kNoisy %s, kPoorCalib %s, kSaturated %s, kLeadingEdgeRecovered %s, kNeighboursRecovered %s, kTowerRecovered %s, kDead %s, kKilled %s, kTPSaturated %s,kL1SpikeFlag %s, kWeird %s, kDiWeird %s, kHasSwitchToGain6 %s, kHasSwitchToGain1 %s, kUnknown %s' %(   evt.kGood[ixtal],  evt.kPoorReco[ixtal],          evt.kOutOfTimE[ixtal],           evt.kFaultyHardware[ixtal],        evt.kNoisy[ixtal],                 evt.kPoorCalib[ixtal],            evt.kSaturated[ixtal],             evt.kLeadingEdgeRecovered[ixtal],  evt.kNeighboursRecovered[ixtal],   evt.kTowerRecovered[ixtal],        evt.kDead[ixtal],                 evt.kKilled[ixtal],                evt.kTPSaturated[ixtal],           evt.kL1SpikeFlag[ixtal],           evt.kWeird[ixtal],                 evt.kDiWeird[ixtal],               evt.kHasSwitchToGain6[ixtal],      evt.kHasSwitchToGain1[ixtal],      evt.kUnknown[ixtal] )
            print evt.runNumber, evt.lumiBlock, evt.eventNumber
        
        



outfile=ROOT.TFile.Open("output_revoXtal_PierreTag_sum8gt0.root", "RECREATE")
hRecovXtalSize.Write()
hGood_recov1.Write()                 
hPoorReco_recov1.Write()             
hOutOfTimE_recov1.Write()            
hFaultyHardware_recov1.Write()       
hNoisy_recov1.Write()                
hPoorCalib_recov1.Write()            
hSaturated_recov1.Write()            
hLeadingEdgeRecovered_recov1.Write() 
hNeighboursRecovered_recov1.Write()  
hTowerRecovered_recov1.Write()       
hDead_recov1.Write()                 
hKilled_recov1.Write()               
hTPSaturated_recov1.Write()          
hL1SpikeFlag_recov1.Write()          
hWeird_recov1.Write()                
hDiWeird_recov1.Write()              
hHasSwitchToGain6_recov1.Write()     
hHasSwitchToGain1_recov1.Write()

outfile.Close()     
