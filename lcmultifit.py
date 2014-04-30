
from numpy import *
from lcspot import *
from lclmleastsq import vincfit
from lclmleastsq import fincfit
from lclmleastsq import paramdist
from lclmleastsq import genrandlc
from lclmleastsq import genfinclc
import sys
from time import sleep
if '-p' in sys.argv: import matplotlib.pyplot as plt

def genmultilc(phase=arange(0,1,.0125), noisefactor=0.02, nspots=2):
  x,intensity,tparams = genrandlc(phase, 0.0)
  tinc = tparams[0]
  tspots = [tparams[1:]]
  for i in range(1,nspots):
    x,lc,tp = genfinclc(tinc, phase, 0.0)
    intensity = lccomb(intensity, lc)
    tspots.append(tp[1:])
  if noisefactor > 0.0:
    lcrange = max(intensity) - min(intensity)
    noise = random.normal(0,lcrange*noisefactor,len(intensity))
    intensity = intensity + noise
  return phase, intensity, tinc, tspots

def sequentialfit(lctuple, nspots, plsprint='some', plsplot=False):
  if nspots != int(nspots) or nspots < 1:
    raise ValueError('nspots must be positive integer')
  
  phase, intensity, tinc, tspots = lctuple
  
  if plsprint != 'none':
    print 'tinc:', tinc, 'tspots:', tspots
    print
  
  fit1ps = vincfit((phase,intensity), plsprint=plsprint, plsplot=plsplot)
  #incs = [p[0] for p in fit1ps]
  #spots = [p[1:] for p in fit1ps]
  
  paramsets = []
  if nspots == 1:
    paramsets = [[p[0],p[1:]] for p in fit1ps]
  else:
    for p in fit1ps:
      inc = p[0]
      firstspot = p[1:]
      iresidual = lcsep(intensity, lcspotfi(phase,inc,firstspot))
      rspots = sequentialfithelper((phase,iresidual),inc,nspots,2,plsprint,plsplot)
      paramsets = paramsets + [[inc,firstspot] + rs for rs in rspots]
  
  fitlcs = [lcmultispot(phase,pset) for pset in paramsets]
  sses = [sum((fitlc - intensity)**2) for fitlc in fitlcs]
  stups = sorted(zip(sses,paramsets,fitlcs), key=lambda x: x[0])
  sses,paramsets,fitlcs = zip(*stups)
  
  if plsprint != 'none':
    print 'tinc:', tinc, 'tspots:', tspots
    print
    print 'best final fit:'
    print paramsets[0]
    print 'sse:', sses[0]
    print 'intensity:'
    print fitlcs[0]
    print
  
  if plsprint == 'all':
    print 'all final fits:'
    for i in range(0,len(paramsets)):
      print paramsets[i], 'sse:', sses[i]
    print
  
  if plsplot:
    plt.figure()
    plt.plot(phase, intensity, 'r.', phase, fitlcs[0], 'b')
    plt.xlabel('phase')
    plt.ylabel('intensity')
    #plt.show()
  
  return paramsets


def sequentialfithelper(lctuple, inc, nspots, spotnum, plsprint='some', plsplot=False):
  phase, intensity = lctuple
  spotps = fincfit(lctuple, inc, plsprint=plsprint, plsplot=plsplot, snum=spotnum)
  if spotnum == nspots:
    return [[p] for p in spotps]
  else:
    spotsets = []
    for p in spotps:
      iresidual = lcsep(intensity, lcspotfi(phase,inc,p))
      sss = sequentialfithelper((phase,iresidual),inc,nspots,spotnum+1,plsprint,plsplot)
      spotsets = spotsets + [[p] + ss for ss in sss]
    return spotsets


def main():
  if '-p' in sys.argv: plt.ion()
  if '-a' in sys.argv:
    pr = 'all'
  else:
    pr = 'some'
  ns = int(sys.argv[-1])
  sequentialfit(genmultilc(nspots=ns),ns,plsprint=pr,plsplot='-p' in sys.argv)
  if '-p' in sys.argv:
    for i in plt.get_fignums():
      plt.figure(i)
      plt.savefig('figure%d.png' % i)
  
if 'lcmultifit.py' in sys.argv:
  main()
