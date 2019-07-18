
runs=[]
lines=[]
lumis=[]
evts=[]


for line in open('eventList_sum8gt0.txt', 'r'):
    myline = line.rstrip('\n')
    lines.append(myline.replace(' ',':'))
    run,lumi,evt = myline.split(' ')
    runs.append(int(run))
    lumis.append(int(lumi))
    evts.append(int(evt))
lines=list(set(lines))
outfile=open('cleanedEventList.txt','w')
for line in lines:
    outfile.write(line+'\n')
    
outfile.close()



runs=list(set(runs))
print runs
lumis=list(set(lumis))
