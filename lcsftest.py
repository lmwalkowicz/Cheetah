
from numpy import *
from lcspot import *
from lcsinglefit import *
from lcgenerate import *
import sys

def main():
  ntrials = int(sys.argv[1])
  isteps = int(sys.argv[2])
  nclusters = int(sys.argv[3])
  threshratio = float(sys.argv[4])
  noisefac = float(sys.argv[5])
  
  for i in range(ntrials):
    phase, intesity, tparams = genrandlc(noisefactor=noisefac)
    paramsets = vincfit((phase,intesity),[],isteps,nclusters,threshratio,plsprint='none',plsplot=False)
    pdists = [paramdist(fps,tparams) for fps in paramsets]
    print pdists

if 'lcsftest.py' in sys.argv:
  main()
