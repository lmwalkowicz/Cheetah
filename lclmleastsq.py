
from numpy import *
from scipy.optimize import leastsq
from lcspot import *
from scipy.cluster.vq import whiten
from scipy.cluster.vq import kmeans2
import sys
if '-p' in sys.argv: import matplotlib.pyplot as plt

def genrandlc(phase=arange(0,1,.025), noisefactor=0.02):
  while True:
    tinc = random.uniform(0,180)
    tlon = random.uniform(0,360)
    tlat = random.uniform(0,90)
    trad = random.uniform(10,20)
    tparams = array([tinc, tlon, tlat, trad])
    intensity = lcspot(phase, tparams)
    if min(intensity) < .999: break
  if noisefactor > 0.0:
    lcrange = max(intensity) - min(intensity)
    noise = random.normal(0,lcrange*noisefactor,len(intensity))
    intensity = intensity + noise
  return phase, intensity, tparams


def genfinclc(tinc, phase=arange(0,1,.025), noisefactor=0.02):
  while True:
    tlon = random.uniform(0,360)
    tlat = random.uniform(0,90)
    trad = random.uniform(10,20)
    tparams = array([tinc, tlon, tlat, trad])
    intensity = lcspot(phase, tparams)
    if min(intensity) < .999: break
  if noisefactor > 0.0:
    lcrange = max(intensity) - min(intensity)
    noise = random.normal(0,lcrange*noisefactor,len(intensity))
    intensity = intensity + noise
  return phase, intensity, tparams


def singlefit(lctuple, plsprint='some', plsplot=False):
  phase, intensity, tparams = lctuple
  
  if plsprint != 'none':
    print 'true params:', tparams
    print 'lightcurve:', intensity
  
  #fit using leastsq, starting in the center of param space 
  fitparams, iem = lcfit(phase,intensity)
  
  if plsprint != 'none':
    print 'fit params:', fitparams
    print 'iem:', iem
  
  #plot
  if plsplot:
    if iem in [1,2,3,4]:
      plt.figure()
      plt.plot(phase, intensity, 'r.', phase, lcspot(phase,fitparams), 'b')
      plt.xlabel('phase')
      plt.ylabel('intensity')
      #plt.show()
    else:
      plt.figure()
      plt.plot(phase, intensity, 'r.')
      plt.xlabel('phase')
      plt.ylabel('intensity')
      #plt.show()
  
  return fitparams


def multifit(lctuple, plsprint='some', plsplot=False):
  phase, intensity, tparams = lctuple
  
  if plsprint != 'none':
    print 'true params:', tparams
    print 'lightcurve:'
    print intensity
    print
  
  #fit using leastsq, starting at various points in param space
  opts,sses,p0s,pdists = lcfits(phase,intensity,tps=tparams)
  
  #if no fits found, plot original curve and exit
  if not opts:
    print 'no fits found =('
    if plsplot:
      plt.figure()
      plt.plot(phase, intensity, 'r.')
      plt.xlabel('phase')
      plt.ylabel('intensity')
      #plt.show()
    return opts, sses, p0s, pdists
  
  #sort by SSE
  sses,opts,p0s,pdists = zip(*sorted(zip(sses,opts,p0s,pdists), key=lambda x: x[0]))
  
  #find best fit
  bestsse, bidx = min((val, idx) for (idx, val) in enumerate(sses))
  bestparams = opts[bidx]
  bestp0s = p0s[bidx]
  bestpdist = pdists[bidx]
  bestlc = lcspot(phase,bestparams)
  
  if plsprint != 'none':
    print 'best fit:'
    print 'opt:',bestparams,'sse:',bestsse,'p0:',bestp0s,'pdist:',bestpdist
    print 'best lightcurve:'
    print bestlc
    print
  
  #print all fits
  if plsprint == 'all':
    print 'all fits:'
    for i in range(0,len(opts)):
      print 'opt:',opts[i],'sse:',sses[i],'p0:',p0s[i],'pdist:',pdists[i]
    print
  
  #plot
  if plsplot:
    plt.figure()
    plt.plot(phase, intensity, 'r.', phase, bestlc, 'b')
    plt.xlabel('phase')
    plt.ylabel('intensity')
    #plt.show()
  
  return bestparams


def smartfit1(lctuple, initmaxfev=20, plsprint='some', plsplot=False):
  phase, intensity, tparams = lctuple
  if plsprint != 'none':
    print 'true params:', tparams
    print 'lightcurve:'
    print intensity
    print
  
  #initial fits using leastsq, starting at various points in param space
  opts,sses,p0s,pdists = lcfits(phase,intensity,tps=tparams,mfev=initmaxfev)
  
  #sort by sses, take best fits
  stups = sorted(zip(sses,opts,p0s,pdists), key=lambda x: x[0])
  sses,opts,p0s,pdists = zip(*stups)
  opts = list(opts)
  
  #print all fits
  #if plsprint == 'all':
  #  print 'all fits:'
  #  for i in range(0,len(opts)):
  #    print 'opt:',opts[i],'sse:',sses[i],'p0:',p0s[i],'pdist:',pdists[i]
  #  print
  
  #full fits
  opts,sses,p1s,pdists = lcfits(phase,intensity,p0stotry=opts[:20],tps=tparams)
  
  #find best fit
  bestsse, bidx = min((val, idx) for (idx, val) in enumerate(sses))
  bestparams = opts[bidx]
  bestp0 = p0s[bidx]
  bestp1 = p1s[bidx]
  bestpdist = pdists[bidx]
  bestlc = lcspot(phase,bestparams)
  
  #print best fit
  if plsprint != 'none':
    print 'best fit:'
    print 'opt:',bestparams,'sse:',bestsse,'p0:',bestp0,'p1:',bestp1,'pdist:',bestpdist
    print 'best lightcurve:'
    print bestlc
    print
  
  #print all fits
  if plsprint == 'all':
    print 'all fits:'
    for i in range(0,len(opts)):
      print 'opt:',opts[i],'sse:',sses[i],'p0:',p0s[i],'p1:',p1s[i],'pdist:',pdists[i]
    print
  
  #plot
  if plsplot:
    plt.figure()
    plt.plot(phase, intensity, 'r.', phase, bestlc, 'b')
    plt.xlabel('phase')
    plt.ylabel('intensity')
    #plt.show()
  
  return bestparams


def smartfit2(lctuple, initmaxfev=20, nclusters=20, plsprint='some', plsplot=False):
  phase, intensity, tparams = lctuple
  
  if plsprint != 'none':
    print 'true params:', tparams
    print 'lightcurve:'
    print intensity
    print
  
  #initial fits using leastsq, starting at various points in param space
  opts,sses,p0s,pdists = lcfits(phase,intensity,tps=tparams,mfev=initmaxfev)
  
  #sort by sses
  stups = sorted(zip(sses,opts,p0s,pdists), key=lambda x: x[0])
  sses,opts,p0s,pdists = zip(*stups)
  
  #cluster best fits
  optsmat = whiten(array(opts[0:3*len(opts)/4]))
  centroid,label = kmeans2(optsmat,nclusters,iter=20,minit='points')
  ll = list(label)
  newp0s = [opts[ll.index(i)] for i in range(0,nclusters) if i in ll]
  
  #full fits
  opts,sses,p1s,pdists = lcfits(phase,intensity,p0stotry=newp0s,tps=tparams)
  
  #sort by sses again
  stups = sorted(zip(sses,opts,p1s,pdists), key=lambda x: x[0])
  sses,opts,p1s,pdists = zip(*stups)
  
  #find best fit(s)
  threshold = sses[0]*2
  bestps = [opts[0]]
  bestsses = [sses[0]]
  bestp1s = [p1s[0]]
  bestpdists = [pdists[0]]
  for i in range(1,len(opts)):
    if sses[i] > threshold:
      break
    if min([paramdist(opts[i],k) for k in bestps]) > 0.001:
      bestps.append(opts[i])
      bestsses.append(sses[i])
      bestp1s.append(p1s[i])
      bestpdists.append(pdists[i])
  
  #print best fit(s)
  if plsprint != 'none':
    print 'best fits:'
    for i in range(0,len(bestps)):
      print 'opt:',bestps[i],'sse:',bestsses[i],'p1:',bestp1s[i],'pdist:',bestpdists[i]
    print
  
  #print all fits
  if plsprint == 'all':
    print 'all fits:'
    for i in range(0,len(opts)):
      print 'opt:',opts[i],'sse:',sses[i],'p1:',p1s[i],'pdist:',pdists[i]
    print
  
  #plot
  if plsplot:
    colors = 'bgcmyk'
    plt.figure()
    plt.plot(phase, intensity, 'r.')
    for i in range(0,len(bestps)):
      plt.plot(phase, lcspot(phase,bestps[i]), colors[i%len(colors)])
    plt.xlabel('phase')
    plt.ylabel('intensity')
    #plt.show()
  
  return bestps


def vincfit(lctuple, initmaxfev=20, nclusters=20, plsprint='some', plsplot=False):
  phase, intensity = lctuple
  
  if plsprint != 'none':
    print 'first spot'
    print 'lightcurve:'
    print intensity
    print
  
  #initial fits using leastsq, starting at various points in param space
  opts,sses,p0s = lcfits(phase,intensity,mfev=initmaxfev)
  
  #sort by sses
  stups = sorted(zip(sses,opts,p0s), key=lambda x: x[0])
  sses,opts,p0s = zip(*stups)
  
  #cluster best fits
  optsmat = whiten(array(opts[0:3*len(opts)/4]))
  centroid,label = kmeans2(optsmat,nclusters,iter=20,minit='points')
  ll = list(label)
  newp0s = [opts[ll.index(i)] for i in range(0,nclusters) if i in ll]
  
  #full fits
  opts,sses,p1s = lcfits(phase,intensity,p0stotry=newp0s)
  
  #sort by sses again
  stups = sorted(zip(sses,opts,p1s), key=lambda x: x[0])
  sses,opts,p1s = zip(*stups)
  
  #find best fit(s)
  threshold = sses[0]*2
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
    plt.xlabel('phase')
    plt.ylabel('intensity')
    #plt.show()
  
  return bestps


def fincfit(lctuple, inc, initmaxfev=20, nclusters=10, plsprint='some', plsplot=False, snum=-1):
  phase, intensity = lctuple
  
  if plsprint != 'none':
    print 'inc:', inc, 'spot #', snum
    print 'residual lightcurve:'
    print intensity
    print
  
  #initial fits using leastsq, starting at various points in param space
  opts,sses,p0s = filcfits(phase,intensity,inc,mfev=initmaxfev)
  
  #sort by sses
  stups = sorted(zip(sses,opts,p0s), key=lambda x: x[0])
  sses,opts,p0s = zip(*stups)
  
  #cluster best fits
  optsmat = whiten(array(opts[0:3*len(opts)/4]))
  centroid,label = kmeans2(optsmat,nclusters,iter=20,minit='points')
  ll = list(label)
  newp0s = [opts[ll.index(i)] for i in range(0,nclusters) if i in ll]
  
  #full fits
  opts,sses,p1s = filcfits(phase,intensity,inc,p0stotry=newp0s)
  
  #sort by sses again
  stups = sorted(zip(sses,opts,p1s), key=lambda x: x[0])
  sses,opts,p1s = zip(*stups)
  
  #find best fit(s)
  threshold = sses[0]*2
  bestps = [opts[0]]
  bestsses = [sses[0]]
  bestp1s = [p1s[0]]
  for i in range(1,len(opts)):
    if sses[i] > threshold:
      break
    if min([paramdist(opts[i],k,array([180.0, 45.0, 7.5])) for k in bestps]) > 0.001:
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


def lcfit(phase,intensity):
  return leastsq(lcspotdiffs3, [90.0,180.0,45.0,7.5], args=(phase,intensity))


def lcfits(phase,intensity,p0stotry=[],tps=[],mfev=0):
  if not p0stotry:
   #p0stotry = [list(x) for x in list(product(range(15,180,30),[180.0],range(7,90,15),range(6,10,1)))]
   for inc in range(15,180,30):
    for lon in [180.0]:
     for lat in range(7,90,15):
      for rad in range(10,20,1):
       p0stotry.append([inc,lon,lat,rad])
  opts = []
  sses = []
  p0s = []
  pdists = []
  for p0 in p0stotry:
    try:
      fps, iem = leastsq(lcspotdiffs3, p0, args=(phase,intensity), maxfev=mfev)
      #fps[1] = fps[1] % 360.0
      fps[3] = abs(fps[3])
      opts.append(fps)
      sses.append(lcspotsse(fps,phase,intensity))
      p0s.append(p0)
      if len(tps) > 0: pdists.append(paramdist(fps,tps))
    except RuntimeError as e:
      print "Runtime Error({0}): {1}".format(e.errno, e.strerror)
  if len(tps) > 0:
    return opts,sses,p0s,pdists
  else:
    return opts,sses,p0s


def filcfits(phase,intensity,inc,p0stotry=[],tps=[],mfev=0):
  if not p0stotry:
   #p0stotry = [list(x) for x in list(product(range(15,180,30),[180.0],range(7,90,15),range(6,10,1)))]
   for lon in [180.0]:
    for lat in range(7,90,15):
     for rad in range(10,20,1):
      p0stotry.append([lon,lat,rad])
  opts = []
  sses = []
  p0s = []
  for p0 in p0stotry:
    try:
      fps, iem = leastsq(lcspotdiffsfi, p0, args=(phase,intensity,inc), maxfev=mfev)
      #fps[0] = fps[1] % 360.0
      fps[2] = abs(fps[2])
      opts.append(fps)
      sses.append(lcspotssefi(fps,phase,intensity,inc))
      p0s.append(p0)
    except RuntimeError as e:
      print "Runtime Error({0}): {1}".format(e.errno, e.strerror)
  return opts,sses,p0s


def paramdistold(fps, tps):
  return sum(((fps/tps)-1)**2)


def paramdist(fps, tps, scalevals=array([90.0, 180.0, 45.0, 7.5])):
  return sum(((fps - tps)/scalevals)**2)


def main():
  if '-nn' in sys.argv:
    nf = 0.0
  else:
    nf = 0.02
  if '-a' in sys.argv:
    pr = 'all'
  else:
    pr = 'some'
  
  if '-1' in sys.argv:
    singlefit(genrandlc(noisefactor=nf), plsplot='-p' in sys.argv)
  elif '-m' in sys.argv:
    multifit(genrandlc(noisefactor=nf), plsprint=pr, plsplot='-p' in sys.argv)
  elif '-s1' in sys.argv:
    smartfit1(genrandlc(noisefactor=nf), plsprint=pr, plsplot='-p' in sys.argv)
  else:
    smartfit2(genrandlc(noisefactor=nf), plsprint=pr, plsplot='-p' in sys.argv)

if 'lclmleastsq.py' in sys.argv:
  main()
