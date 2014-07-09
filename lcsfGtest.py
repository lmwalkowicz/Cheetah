
from numpy import *
from lcspot import *
from lcgenerate import *
import lcsinglefit
import lcsinglefitG
import sys

def main():
  isteps = int(sys.argv[1])         # def 20-30
  nclusters = int(sys.argv[2])      # def 20-30
  threshratio = float(sys.argv[3])  # def 2.0
  noisefac = float(sys.argv[4])     # def .002
  meths = ('Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG', 'Anneal')
  
  phase, intesity, tparams = genrandlc(noisefactor=noisefac)
  
  print 'true parameters', tparams
  print
  
  print 'L-M :'
  print
  paramsets = lcsinglefit.vincfit((phase,intesity),[],isteps,nclusters,threshratio,plsprint='all',plsplot=False)
  pdists = [paramdist(fps,tparams) for fps in paramsets]
  print pdists
  print
  
  for meth in meths:
    print meth, ':'
    print
    paramsets = lcsinglefit2.vincfit((phase,intesity),[],isteps,nclusters,threshratio,meth,plsprint='all',plsplot=False)
    pdists = [paramdist(fps,tparams) for fps in paramsets]
    print pdists
    print

if 'lcsf2test.py' in sys.argv:
  main()
