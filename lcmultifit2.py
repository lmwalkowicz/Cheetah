
from numpy import *
from lcspot2 import *
from lcsinglefit2 import *
import matplotlib.pyplot as plt
plt.ion()


def doublefit(lctuple, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some', plsplot=False):
  time, intensity = lctuple
  fit1ps = vincfit(lctuple, p0s, initsteps, nclusters, threshratio, plsprint=plsprint, plsplot=plsplot)
  paramsets = []
  for p1 in fit1ps:
    inc = p1[0]
    iresidual = lcsep(intensity, lcspot(time,p1))
    fit2ps = fincfit((time,iresidual), inc, p0s, initsteps, nclusters, threshratio, plsprint=plsprint, plsplot=plsplot)
    for p2 in fit2ps:
      lat1 = p1[2]
      T1 = p1[4]
      lat2 = p2[1]
      T2 = p2[3]
      alpha = (T1 - T2)/((T1 * sin(lat1)**2) - (T2 * sin(lat2)**2))
      Teq = T1 * (1 - alpha * sin(lat1)**2)
      paramsets.append([inc, Teq, alpha, p1[1:4], p2[0:3]])
  return simulfit(lctuple, 2, paramsets, threshratio, plsprint, plsplot)


def simulfit(lctuple, nspots, p0s, threshratio=2, plsprint='some', plsplot=False):
  if nspots != int(nspots) or nspots < 1:
    raise ValueError('nspots must be positive integer')
  time, intensity = lctuple
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
      fps, iem = leastsq(lcmultispotdiffs, p0, args=(time,intensity))
      opts.append(fps)
      sses.append(lcmultispotsse(fps,time,intensity))
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
    plt.xlabel('time')
    plt.ylabel('intensity')
    plt.title('simultaneous ' + str(nspots) + ' fits')
    #plt.show()
  
  return bestps


def ratchetfit(lctuple, nspots, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some', plsplot=False):
  if nspots != int(nspots) or nspots < 2:
    raise ValueError('nspots must be positive integer >= 2')
  time, intensity = lctuple
  paramsets = doublefit(lctuple, p0s, initsteps, nclusters, threshratio, plsprint, plsplot)
  for spotnum in range(3,nspots+1):
    newparamsets = []
    for pset in paramsets:
      iresidual = lcsep(intensity, lcmultispot(time,pset))
      spotps = fstarfit((time,iresidual), pset[0], pset[1], pset[2], [], initsteps, nclusters, threshratio, plsprint, plsplot, spotnum)
      newparamsets = newparamsets + [pset + [spotp] for spotp in spotps]
    paramsets = simulfit(lctuple, spotnum, newparamsets, threshratio, plsprint, plsplot)
  return paramsets

