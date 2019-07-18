import ROOT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='infilename', action='store',help='input file.root', required=True)
parser.add_argument('--output', dest='outtxt', action='store', help='output file.txt', required=True)
args = parser.parse_args()

infile=ROOT.TFile.Open(args.infilename, "READ")

ntuple=infile.Get("ntupler/selected")

outfile=open(args.outtxt, "w")

eventlist=[]
for row in ntuple:
    run=row.runNumber
    lumi=row.lumiBlock
    evt=row.eventNumber
    eventlist.append("%d:%d:%d" %(run, lumi, evt))

eventlist=list(set(eventlist))

for line in eventlist :
    outfile.write(line+"\n")

#outfile.write("%d:%d:%d\n" %(run, lumi, evt))
    
outfile.close()


    

