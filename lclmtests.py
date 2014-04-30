
from numpy import *
from lcspot import *
from lclmleastsq import *
import sys
if '-p' in sys.argv: import matplotlib.pyplot as plt

def testfits(fitalg, nf=0.02, iters=100):
  pdists = []
  for r in range(iters):
    rlc = genrandlc(noisefactor=nf)
    tps = rlc[2]
    fps = fitalg(rlc, plsprint='none', plsplot=False)
    pdists.append(paramdist(fps,tps))
  return pdists

def main():
  if '-nn' in sys.argv:
    nf = 0.0
  else:
    nf = 0.02
  
  if '-1' in sys.argv:
    fitalg = singlefit
    name = 'singlefit'
  elif '-m' in sys.argv:
    fitalg = multifit
    name = 'multifit'
  elif '-s1' in sys.argv:
    fitalg = smartfit1
    name = 'smartfit1'
  elif '-s2' in sys.argv:
    fitalg = smartfit2
    name = 'smartfit2'
  else:
    return
  
  pdists = testfits(fitalg, nf, 5)
  
  print 'pdists:', pdists
  print 'mean:', mean(pdist)
  
  if '-p' in sys.argv:
    plt.hist(pdists)
    plt.title(name)
    plt.xlabel('Scaled Parameter Distance')
    plt.ylabel('Frequency')
    plt.show()

if 'lclmtests.py' in sys.argv:
  main()
