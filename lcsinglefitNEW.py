
from numpy import *
from scipy.optimize import leastsq
from lcspotNEW import *
from scipy.cluster.vq import whiten
from scipy.cluster.vq import kmeans2

def vincfit(lctuple, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some'):
  time, intensity = lctuple
  
  if plsprint != 'none':
    print 'vincfit with lightcurve:'
    print intensity
    print 'p0s = '
    print p0s
    print
  
  #initial fits using leastsq, starting at various points in param space
  opts,sses,p0s = vilcfits(time,intensity,p0s,initsteps)
  opts = [boundparams(p) for p in opts]
  
  #sort by sses
  stups = sorted(zip(sses,opts,p0s), key=lambda x: x[0])
  sses,opts,p0s = zip(*stups)
  
  #cluster best fits
  optsmat = whiten(array(opts[0:3*len(opts)/4]))
  centroid,label = kmeans2(optsmat,nclusters,iter=20,minit='points')
  ll = list(label)
  newp0s = [opts[ll.index(i)] for i in range(0,nclusters) if i in ll]
  
  #full fits
  opts,sses,p1s = vilcfits(time,intensity,newp0s)
  opts = [boundparams(p) for p in opts]
  
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
    if min([paramdist(opts[i],k,scalevals=array([45.0, 180.0, 90.0, 17.5, 25.0])) for k in bestps]) > 0.001:
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
  
  return bestps


def fincfit(lctuple, inc, p0s=[], initsteps=20, nclusters=20, threshratio=2, plsprint='some'):
  time, intensity = lctuple
  
  if plsprint != 'none':
    print 'fixed-inc fit, inc:', inc
    print 'residual lightcurve:'
    print intensity
    print 'p0s = '
    print p0s
    print
  
  #initial fits using leastsq, starting at various points in param space
  opts,sses,p0s = filcfits(time,intensity,inc,p0s,initsteps)
  opts = [boundparamsfi(p) for p in opts]
  
  #sort by sses
  stups = sorted(zip(sses,opts,p0s), key=lambda x: x[0])
  sses,opts,p0s = zip(*stups)
  
  #cluster best fits
  optsmat = whiten(array(opts[0:3*len(opts)/4]))
  centroid,label = kmeans2(optsmat,nclusters,iter=20,minit='points')
  ll = list(label)
  newp0s = [opts[ll.index(i)] for i in range(0,nclusters) if i in ll]
  
  #full fits
  opts,sses,p1s = filcfits(time,intensity,inc,newp0s)
  opts = [boundparamsfi(p) for p in opts]
  
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
    if min([paramdist(opts[i],k,scalevals=array([180.0, 90.0, 17.5, 25.0])) for k in bestps]) > 0.001:
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
  
  return bestps


def fstarfit(lctuple, inc, teq, alpha, p0s=[], initsteps=20, nclusters=20, threshratio=2, plsprint='some', snum=-1):
  time, intensity = lctuple
  
  if plsprint != 'none':
    print 'fixed-star fit, inc:', inc, 'teq:', teq, 'alpha:', alpha, 'spot #', snum
    print 'residual lightcurve:'
    print intensity
    print 'p0s = '
    print p0s
    print
  
  #initial fits using leastsq, starting at various points in param space
  opts,sses,p0s = fstarfits(time,intensity,inc,teq,alpha,p0s,initsteps)
  opts = [boundparamsfstar(p) for p in opts]
  
  #sort by sses
  stups = sorted(zip(sses,opts,p0s), key=lambda x: x[0])
  sses,opts,p0s = zip(*stups)
  
  #cluster best fits
  optsmat = whiten(array(opts[0:3*len(opts)/4]))
  centroid,label = kmeans2(optsmat,nclusters,iter=20,minit='points')
  ll = list(label)
  newp0s = [opts[ll.index(i)] for i in range(0,nclusters) if i in ll]
  
  #full fits
  opts,sses,p1s = fstarfits(time,intensity,inc,teq,alpha,newp0s)
  opts = [boundparamsfstar(p) for p in opts]
  
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
    if min([paramdist(opts[i],k,scalevals=array([180.0, 90.0, 17.5])) for k in bestps]) > 0.001:
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
  
  return bestps


def vilcfits(time,intensity,p0s=[],mfev=0):
  if not p0s:
    p0s = spacedparams(5)
  opts = []
  sses = []
  p0s2 = []
  for p0 in p0s:
    try:
      fps, iem = leastsq(lcspotdiffs3, p0, args=(time,intensity), maxfev=mfev)
      opts.append(fps)
      sses.append(lcspotsse3(fps,time,intensity))
      p0s2.append(p0)
    except RuntimeError as e:
      print "Runtime Error({0}): {1}".format(e.errno, e.strerror)
  return opts,sses,p0s2


def filcfits(time,intensity,inc,p0s=[],mfev=0):
  if not p0s:
    p0s = spacedparamsfi(5)
  opts = []
  sses = []
  p0s2 = []
  for p0 in p0s:
    try:
      fps, iem = leastsq(lcspotdiffsfi3, p0, args=(time,intensity,inc), maxfev=mfev)
      opts.append(fps)
      sses.append(lcspotssefi(fps,time,intensity,inc))
      p0s2.append(p0)
    except RuntimeError as e:
      print "Runtime Error({0}): {1}".format(e.errno, e.strerror)
  return opts,sses,p0s2


def fstarfits(time,intensity,inc,teq,alpha,p0s=[],mfev=0):
  if not p0s:
    p0s = spacedparamsspot(5)
  opts = []
  sses = []
  p0s2 = []
  for p0 in p0s:
    try:
      fps, iem = leastsq(lcspotdiffsfstar3, p0, args=(time,intensity,inc,teq,alpha), maxfev=mfev)
      opts.append(fps)
      sses.append(lcspotssefstar(fps,time,intensity,inc,teq,alpha))
      p0s2.append(p0)
    except RuntimeError as e:
      print "Runtime Error({0}): {1}".format(e.errno, e.strerror)
  return opts,sses,p0s2
