
from numpy import *
from lcspot import *
from lcsinglefit import *
from lcgenerate import *
import sys
if '-p' in sys.argv: import matplotlib.pyplot as plt


def sequentialfit(lctuple, nspots, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some', plsplot=False):
  if nspots != int(nspots) or nspots < 1:
    raise ValueError('nspots must be positive integer')
  
  phase, intensity = lctuple
  
  fit1ps = vincfit((phase,intensity), p0s, initsteps, nclusters, threshratio, plsprint=plsprint, plsplot=plsplot)
  
  paramsets = []
  if nspots == 1:
    paramsets = [[p[0],p[1:]] for p in fit1ps]
  else:
    for p in fit1ps:
      inc = p[0]
      firstspot = p[1:]
      iresidual = lcsep(intensity, lcspotfi(phase,inc,firstspot))
      rspots = sequentialfithelper((phase,iresidual),inc,nspots,2,p0s,initsteps,nclusters,threshratio,plsprint,plsplot)
      paramsets = paramsets + [[inc,firstspot] + rs for rs in rspots]
  
  fitlcs = [lcmultispot(phase,pset) for pset in paramsets]
  sses = [sum((fitlc - intensity)**2) for fitlc in fitlcs]
  stups = sorted(zip(sses,paramsets,fitlcs), key=lambda x: x[0])
  sses,paramsets,fitlcs = zip(*stups)
  
  if plsprint != 'none':
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


def sequentialfithelper(lctuple, inc, nspots, spotnum, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some', plsplot=False):
  phase, intensity = lctuple
  spotps = fincfit(lctuple, inc, p0s, initsteps, nclusters, threshratio, plsprint, plsplot, spotnum)
  if spotnum == nspots:
    return [[p] for p in spotps]
  else:
    spotsets = []
    for p in spotps:
      iresidual = lcsep(intensity, lcspotfi(phase,inc,p))
      sss = sequentialfithelper((phase,iresidual),inc,nspots,spotnum+1,p0s,initsteps,nclusters,threshratio,plsprint,plsplot)
      spotsets = spotsets + [[p] + ss for ss in sss]
    return spotsets


def simulfit(lctuple, nspots, p0s, threshratio=2, plsprint='some', plsplot=False):
  if nspots != int(nspots) or nspots < 1:
    raise ValueError('nspots must be positive integer')
  phase, intensity = lctuple
  if plsprint != 'none':
    print 'lightcurve:'
    print intensity
    print
  
  flatp0s = [gpl2fpl(p0) for p0 in p0s]
  opts = []
  sses = []
  p0s2 = []
  for p0 in flatp0s:
    try:
      fps, iem = leastsq(lcmultispotdiffs, p0, args=(phase,intensity))
      opts.append(fps)
      sses.append(lcmultispotsse(fps,phase,intensity))
      p0s2.append(p0)
    except RuntimeError as e:
      print "Runtime Error({0}): {1}".format(e.errno, e.strerror)
  
  stups = sorted(zip(sses,opts,p0s2), key=lambda x: x[0])
  sses,opts,p0s2 = zip(*stups)
  
  opts = [fpl2gpl(p) for p in opts]
  p0s2 = [fpl2gpl(p) for p in p0s2]
  
  threshold = sses[0]*threshratio
  bestps = [opts[0]]
  bestsses = [sses[0]]
  bestp0s = [p0s2[0]]
  for i in range(1,len(opts)):
    if sses[i] > threshold:
      break
    bestps.append(opts[i])
    bestsses.append(sses[i])
    bestp0s.append(p0s2[i])
  
  #print best fit(s)
  if plsprint != 'none':
    print 'best fits:'
    for i in range(0,len(bestps)):
      print 'opt:',bestps[i],'sse:',bestsses[i],'p0:',bestp0s[i]
    print
  
  #print all fits
  if plsprint == 'all':
    print 'all fits:'
    for i in range(0,len(opts)):
      print 'opt:',opts[i],'sse:',sses[i],'p0:',p0s2[i]
    print
  
  #plot
  if plsplot:
    colors = 'bgcmyk'
    plt.figure()
    plt.plot(phase, intensity, 'r.')
    for i in range(0,len(bestps)):
      plt.plot(phase, lcmultispot(phase,bestps[i]), colors[i%len(colors)])
    plt.xlabel('phase')
    plt.ylabel('residual intensity')
    plt.title('simultaneous ' + str(nspots) + ' fits')
    #plt.show()
  
  return bestps


def combinedfit(lctuple, nspots, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some', plsplot=False):
  if nspots != int(nspots) or nspots < 1:
    raise ValueError('nspots must be positive integer')
  paramsets = sequentialfit(lctuple, nspots, p0s, initsteps, nclusters, threshratio, plsprint, plsplot)
  paramsets = simulfit(lctuple, nspots, paramsets, threshratio, plsprint, plsplot)
  return paramsets


def ratchetfit(lctuple, nspots, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some', plsplot=False):
  if nspots != int(nspots) or nspots < 2:
    raise ValueError('nspots must be positive integer >= 2')
  phase, intensity = lctuple
  paramsets = combinedfit(lctuple, 2, p0s, initsteps, nclusters, threshratio, plsprint, plsplot)
  for spotnum in range(3,nspots+1):
    newparamsets = []
    for pset in paramsets:
      iresidual = lcsep(intensity, lcmultispot(phase,pset))
      spotps = fincfit((phase,iresidual), pset[0], p0s, initsteps, nclusters, threshratio, plsprint, plsplot, spotnum)
      newparamsets = newparamsets + [pset + [spotp] for spotp in spotps]
    paramsets = simulfit(lctuple, spotnum, newparamsets, threshratio, plsprint, plsplot)
  return paramsets


def main():
  if '-a' in sys.argv:
    pr = 'all'
  else:
    pr = 'some'
  ns = int(sys.argv[-1])
  
  phase, intesity, tparams = genmultilc(nspots=ns,noisefactor=0.0)
  print 'tparams:', tparams
  print
  
  paramsets = ratchetfit((phase,intesity),ns,plsprint=pr,plsplot='-p' in sys.argv)
  pdists = [multiparamdists(fps,tparams) for fps in paramsets]
  
  print
  print 'tparams:', tparams
  print 'final param sets:'
  for i in range(len(paramsets)):
    print 'params: ', paramsets[i]
    print 'pdists: ', pdists[i]
  
  if '-p' in sys.argv:
    for i in plt.get_fignums():
      plt.figure(i)
      plt.savefig('figure%d.png' % i)

  
if 'lcmultifit.py' in sys.argv:
  main()
