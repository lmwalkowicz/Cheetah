
from numpy import *
from lcspot import *
from lcmultifit import *
from lcgenerate import *
import sys

def main():
  ntrials = int(sys.argv[1])
  nspots = int(sys.argv[2])
  isteps = int(sys.argv[3])
  nclusters = int(sys.argv[4])
  threshratio = float(sys.argv[5])
  noisefac = float(sys.argv[6])
  
  for i in range(ntrials):
    phase, intesity, tparams = genmultilc(nspots=nspots,noisefactor=noisefac)
    paramsets = ratchetfit((phase,intesity),nspots,[],isteps,nclusters,threshratio,plsprint='none',plsplot=False)
    pdists = [multiparamdists(fps,tparams) for fps in paramsets]
    print pdists

if 'lcmftest.py' in sys.argv:
  main()
