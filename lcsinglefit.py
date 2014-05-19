
from numpy import *
from scipy.optimize import leastsq
from lcspot import *
from scipy.cluster.vq import whiten
from scipy.cluster.vq import kmeans2
import sys
if '-p' in sys.argv: import matplotlib.pyplot as plt


def vincfit(lctuple, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some', plsplot=False):
  phase, intensity = lctuple
  
  if plsprint != 'none':
    print 'first spot'
    print 'lightcurve:'
    print intensity
    print
  
  #initial fits using leastsq, starting at various points in param space
  opts,sses,p0s = vilcfits(phase,intensity,p0s,initsteps)
  
  #sort by sses
  stups = sorted(zip(sses,opts,p0s), key=lambda x: x[0])
  sses,opts,p0s = zip(*stups)
  
  #cluster best fits
  optsmat = whiten(array(opts[0:3*len(opts)/4]))
  centroid,label = kmeans2(optsmat,nclusters,iter=20,minit='points')
  ll = list(label)
  newp0s = [opts[ll.index(i)] for i in range(0,nclusters) if i in ll]
  
  #full fits
  opts,sses,p1s = vilcfits(phase,intensity,newp0s)
  
  #sort by sses again
  stups = sorted(zip(sses,opts,p1s), key=lambda x: x[0])
  sses,opts,p1s = zip(*stups)
  
  #find best fit(s)
  threshold = sses[0]*threshratio
  bestps = [opts[0]]
  bestsses = [sses[0]]
  bestp1s = [p1s[0]]
  for i in range(1,len(opts)):
    if sses[i] > threshold:
      break
    if min([paramdist(opts[i],k) for k in bestps]) > 0.001:
      bestps.append(opts[i])
      bestsses.append(sses[i])
      bestp1s.append(p1s[i])
  
  #print best fit(s)
  if plsprint != 'none':
    print 'best fits:'
    for i in range(0,len(bestps)):
      print 'opt:',bestps[i],'sse:',bestsses[i],'p1:',bestp1s[i]
    print
  
  #print all fits
  if plsprint == 'all':
    print 'all fits:'
    for i in range(0,len(opts)):
      print 'opt:',opts[i],'sse:',sses[i],'p1:',p1s[i]
    print
  
  #plot
  if plsplot:
    colors = 'bgcmyk'
    plt.figure()
    plt.plot(phase, intensity, 'r.')
    for i in range(0,len(bestps)):
      plt.plot(phase, lcspot(phase,bestps[i]), colors[i%len(colors)])
    plt.title('initial single spot fit')
    plt.xlabel('phase')
    plt.ylabel('intensity')
    #plt.show()
  
  return bestps


def fincfit(lctuple, inc, p0s=[], initsteps=20, nclusters=20, threshratio=2, plsprint='some', plsplot=False, snum=-1):
  phase, intensity = lctuple
  
  if plsprint != 'none':
    print 'inc:', inc, 'spot #', snum
    print 'residual lightcurve:'
    print intensity
    print
  
  #initial fits using leastsq, starting at various points in param space
  opts,sses,p0s = filcfits(phase,intensity,inc,p0s,initsteps)
  
  #sort by sses
  stups = sorted(zip(sses,opts,p0s), key=lambda x: x[0])
  sses,opts,p0s = zip(*stups)
  
  #cluster best fits
  optsmat = whiten(array(opts[0:3*len(opts)/4]))
  centroid,label = kmeans2(optsmat,nclusters,iter=20,minit='points')
  ll = list(label)
  newp0s = [opts[ll.index(i)] for i in range(0,nclusters) if i in ll]
  
  #full fits
  opts,sses,p1s = filcfits(phase,intensity,inc,newp0s)
  
  #sort by sses again
  stups = sorted(zip(sses,opts,p1s), key=lambda x: x[0])
  sses,opts,p1s = zip(*stups)
  
  #find best fit(s)
  threshold = sses[0]*threshratio
  bestps = [opts[0]]
  bestsses = [sses[0]]
  bestp1s = [p1s[0]]
  for i in range(1,len(opts)):
    if sses[i] > threshold:
      break
    if min([spotparamdist(opts[i],k) for k in bestps]) > 0.001:
      bestps.append(opts[i])
      bestsses.append(sses[i])
      bestp1s.append(p1s[i])
  
  #print best fit(s)
  if plsprint != 'none':
    print 'best fits:'
    for i in range(0,len(bestps)):
      print 'opt:',bestps[i],'sse:',bestsses[i],'p1:',bestp1s[i]
    print
  
  #print all fits
  if plsprint == 'all':
    print 'all fits:'
    for i in range(0,len(opts)):
      print 'opt:',opts[i],'sse:',sses[i],'p1:',p1s[i]
    print
  
  #plot
  if plsplot:
    colors = 'bgcmyk'
    plt.figure()
    plt.plot(phase, intensity, 'r.')
    for i in range(0,len(bestps)):
      plt.plot(phase, lcspotfi(phase,inc,bestps[i]), colors[i%len(colors)])
    plt.xlabel('phase')
    plt.ylabel('residual intensity')
    plt.title('inc: ' + str(inc) + ', spot #' + str(snum))
    #plt.show()
  
  return bestps


def vilcfits(phase,intensity,p0s=[],mfev=0):
  if not p0s:
    p0s = spacedparams(5)
  opts = []
  sses = []
  p0s2 = []
  for p0 in p0s:
    try:
      fps, iem = leastsq(lcspotdiffs3, p0, args=(phase,intensity), maxfev=mfev)
      #fps[1] = fps[1] % 360.0
      #fps[3] = abs(fps[3])
      opts.append(fps)
      sses.append(lcspotsse(fps,phase,intensity))
      p0s2.append(p0)
    except RuntimeError as e:
      print "Runtime Error({0}): {1}".format(e.errno, e.strerror)
  return opts,sses,p0s2


def filcfits(phase,intensity,inc,p0s=[],mfev=0):
  if not p0s:
    p0s = spacedparamsfi(5)
  opts = []
  sses = []
  p0s2 = []
  for p0 in p0s:
    try:
      fps, iem = leastsq(lcspotdiffsfi, p0, args=(phase,intensity,inc), maxfev=mfev)
      #fps[0] = fps[1] % 360.0
      #fps[2] = abs(fps[2])
      opts.append(fps)
      sses.append(lcspotssefi(fps,phase,intensity,inc))
      p0s2.append(p0)
    except RuntimeError as e:
      print "Runtime Error({0}): {1}".format(e.errno, e.strerror)
  return opts,sses,p0s2

