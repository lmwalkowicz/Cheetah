
from numpy import *
from lcspot3 import *
from lcsinglefit3 import *
import matplotlib.pyplot as plt
plt.ion()


def doublefit1(lctuple, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some', plsplot=False):
  time, intensity = lctuple
  fit1ps = vincfit(lctuple, p0s, initsteps, nclusters, threshratio, plsprint=plsprint, plsplot=plsplot)
  paramsets = []
  for p1 in fit1ps:
    inc = p1[0]
    iresidual = lcsep(intensity, lcspot(time,p1))
    iresidual = iresidual/(sorted(iresidual)[int(len(iresidual)*.98)]) #renormalize residual
    fit2ps = fincfit((time,iresidual), inc, p0s, initsteps, nclusters, threshratio, plsprint=plsprint, plsplot=plsplot)
    for p2 in fit2ps:
      lat1 = p1[2]
      T1 = p1[4]
      lat2 = p2[1]
      T2 = p2[3]
      alpha = (T1 - T2)/((T1 * sin(lat1)**2) - (T2 * sin(lat2)**2))
      Teq = T1 * (1 - alpha * sin(lat1)**2)
      paramsets.append([inc, Teq, alpha, p1[1:4], p2[0:3]])
  return paramsets


def doublefit2(lctuple, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some', plsplot=False):
  psets = doublefit1(lctuple, p0s, initsteps, nclusters, threshratio, plsprint, plsplot)
  return simulfit(lctuple, 2, psets, threshratio, plsprint, plsplot)


def simulfit(lctuple, nspots, p0s, threshratio=2, plsprint='some', plsplot=False):
  if nspots != int(nspots) or nspots < 1:
    raise ValueError('nspots must be positive integer')
  time, intensity = lctuple
  if plsprint != 'none':
    print 'lightcurve:'
    print intensity
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
  
  #plot
  if plsplot:
    colors = 'bgcmyk'
    plt.figure()
    plt.plot(time, intensity, 'r.')
    for i in range(0,len(bestps)):
      plt.plot(time, lcmultispot(time,bestps[i]), colors[i%len(colors)])
    plt.xlabel('time')
    plt.ylabel('intensity')
    plt.title('simultaneous ' + str(nspots) + '-spot fits')
    #plt.show()
  
  return bestps


