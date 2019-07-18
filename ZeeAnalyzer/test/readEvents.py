#!/usr/bin/env python

import sys
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()
from DataFormats.FWLite import Events, Handle


events = Events("test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_PU_sum8gt0.root")
txtfile=open("events_PU_sum8gt0.txt","w")
for event_nr,event in enumerate(events):
    runnr = event.eventAuxiliary().run()
    eventnr = event.eventAuxiliary().event()
    lumi = event.eventAuxiliary().luminosityBlock()

    print (runnr,lumi, eventnr)
 
    event_id=str(event.eventAuxiliary().run())+":"+str(event.eventAuxiliary().luminosityBlock())+":"+str(event.eventAuxiliary().event())
    print event_id
    txtfile.write(event_id+"\n")



txtfile.close()
