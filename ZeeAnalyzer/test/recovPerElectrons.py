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

#tree0 =uproot.open("electronTreeZEE_PierreTag_sum8gt20_xtalInfo.root")['ntupler']['selected']
tree0 =uproot.open("electronTreeZEE_PierreTag_sum8gt0_xtalInfo.root")['ntupler']['selected']

myTree0=tree0.arrays(['runNumber', 'lumiBlock', 'eventNumber',  'e1IsRecovered', 'e2IsRecovered', 'e1SeedIEta', 'e2SeedIEta', 'e1SeedIPhi', 'e2SeedIPhi', 'xtalRawId', 'xtalEn', 'vSum8'])


print myTree0['e1IsRecovered']

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

df=pd.DataFrame(np.array([run0, lumi0, evt0, rawId0_0,rawId0_1, rawId0_2, rawId0_3, rawId0_4, myTree0['e1IsRecovered'], myTree0['e2IsRecovered'],  myTree0['e1SeedIEta'], myTree0['e2SeedIEta'], myTree0['e1SeedIPhi'], myTree0['e2SeedIPhi']]).T, columns=['runNumber', 'lumiBlock', 'eventNumber', 'xtalRawId0', 'xtalRawId1', 'xtalRawId2', 'xtalRawId3', 'xtalRawId4',  'e1IsRecovered', 'e2IsRecovered', 'e1SeedIEta', 'e2SeedIEta', 'e1SeedIPhi', 'e2SeedIPhi' ] )

#print df

print df.shape[0], df[df.e1IsRecovered==df.e2IsRecovered].shape[0]
df=df[~(df.e1IsRecovered==df.e2IsRecovered)]

print 'total number of Zee', df.shape[0]
oneRecov=df[ (df.xtalRawId4==0) & (df.xtalRawId3==0) & (df.xtalRawId2==0) & ((df.xtalRawId1==0) | (df.xtalRawId1==df.xtalRawId0))]

print 'is one recov', oneRecov.shape[0]
moreRecov=pd.merge(df, oneRecov, how='outer', on=None, left_on=None, right_on=None,
         left_index=False, right_index=False, sort=True,
         suffixes=('_x', '_y'), copy=True, indicator=True,
         validate=None)
print moreRecov[moreRecov._merge=='left_only'].shape[0]

moreRecov=moreRecov[moreRecov._merge=='left_only']
moreRecov=moreRecov.drop(columns='_merge')
#print moreRecov
def que(x) :
    mylist=[x.xtalRawId4, x.xtalRawId3, x.xtalRawId2, x.xtalRawId1, x.xtalRawId0]
    return len(list(set(mylist)))

moreRecov['nRecov']=moreRecov.apply(que, axis=1)
twoRecov=moreRecov[moreRecov.nRecov==2]
print 'two recov in one electron', twoRecov.shape[0]

threeRecov=moreRecov[moreRecov.nRecov==3]
fourRecov=moreRecov[moreRecov.nRecov==4]
fiveRecov=moreRecov[moreRecov.nRecov==5]
print df.shape[0], oneRecov.shape[0], twoRecov.shape[0], threeRecov.shape[0], fourRecov.shape[0], fiveRecov.shape[0]
