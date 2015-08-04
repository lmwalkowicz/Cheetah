
from numpy import *
from lcspotNEW import *
from lcsinglefitNEW import *

def doublefit1(lctuple, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some'):
  time, intensity = lctuple
  fit1ps = vincfit(lctuple, p0s, initsteps, nclusters, threshratio, plsprint=plsprint)
  paramsets = []
  for p1 in fit1ps:
    inc = p1[0]
    iresidual = lcsep(intensity, lcspot(time,p1))
    iresidual = iresidual/(sorted(iresidual)[int(len(iresidual)*.98)]) #renormalize residual
    fit2ps = fincfit((time,iresidual), inc, p0s, initsteps, nclusters, threshratio, plsprint=plsprint)
    for p2 in fit2ps:
      lat1 = p1[2]
      T1 = p1[4]
      lat2 = p2[1]
      T2 = p2[3]
      alpha = (T1 - T2)/((T1 * sin(lat1)**2) - (T2 * sin(lat2)**2))
      Teq = T1 * (1 - alpha * sin(lat1)**2)
      paramsets.append([inc, Teq, alpha, p1[1:4], p2[0:3]])
  return paramsets


def doublefit2(lctuple, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some'):
  psets = doublefit1(lctuple, p0s, initsteps, nclusters, threshratio, plsprint)
  return simulfit(lctuple, psets, threshratio, plsprint)


def ratchetfit(lctuple, nspots, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some'):
  if nspots != int(nspots) or nspots < 2:
    raise ValueError('nspots must be positive integer >= 2')
  time, intensity = lctuple
  paramsets = doublefit2(lctuple, p0s, initsteps, nclusters, threshratio, plsprint)
  for spotnum in range(3,nspots+1):
    newparamsets = []
    for pset in paramsets:
      iresidual = lcsep(intensity, lcmultispot(time,pset))
      spotps = fstarfit((time,iresidual), pset[0], pset[1], pset[2], [], initsteps, nclusters, threshratio, plsprint, spotnum)
      newparamsets = newparamsets + [pset + [spotp] for spotp in spotps]
    paramsets = simulfit(lctuple, newparamsets, threshratio, plsprint)
  return paramsets


def simulfit(lctuple, p0s, threshratio=2, plsprint='some'):
  time, intensity = lctuple
  if plsprint != 'none':
    print
    print 'simulfit with lightcurve:'
    print intensity
    print 'initial points:'
    for i in range(len(p0s)):
      print p0s[i]
    print
  
  opts = []
  sses = []
  p0s2 = []
  for p0 in p0s:
    try:
      fps, iem = leastsq(lcmultispotdiffs, gpl2fpl(p0), args=(time,intensity))
      opts.append(fps)
      sses.append(lcmultispotsse(fps,time,intensity))
      p0s2.append(p0)
    except RuntimeError as e:
      print "Runtime Error({0}): {1}".format(e.errno, e.strerror)
  
  stups = sorted(zip(sses,opts,p0s2), key=lambda x: x[0])
  sses,opts,p0s2 = zip(*stups)
  
  opts = [boundparamsms(fpl2gpl(p)) for p in opts]
  
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
  
  return bestps

