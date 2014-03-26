
from numpy import *
from scipy.optimize import leastsq
from lcspot import *
import matplotlib.pyplot as plt
import sys

def genrandlc(phase=arange(0,1,.025), noisefactor=0.05):
  while True:
    tinc = random.uniform(0,180)
    tlon = random.uniform(0,360)
    tlat = random.uniform(0,90)
    trad = random.uniform(5,10)
    tparams = array([tinc, tlon, tlat, trad])
    intensity = lcspot(phase, tparams)
    if min(intensity) < .999: break
  if noisefactor > 0.0:
    lcrange = max(intensity) - min(intensity)
    noise = random.normal(0,lcrange*noisefactor,len(intensity))
    intensity = intensity + noise
  return intensity, tparams, phase


def randfit(plsplot=False):
  #generate random single-spot lightcurve
  intensity, tparams, phase = genrandlc()
  print 'true params:', tparams
  print 'lightcurve:', intensity
  
  #fit using leastsq, starting in the center of param space 
  fitparams, iem = lcfit(phase,intensity)
  print 'fit params:', fitparams
  print 'iem:', iem
  
  #plot
  if plsplot:
    if iem in [1,2,3,4]:
      plt.plot(phase, intensity, 'r.', phase, lcspot(phase,fitparams), 'b')
      plt.xlabel('phase')
      plt.ylabel('intensity')
      plt.show()
    else:
      plt.plot(phase, intensity, 'r.')
      plt.xlabel('phase')
      plt.ylabel('intensity')
      plt.show()
  
  return intensity, phase, tparams, fitparams, iem


def randfit2(printall=False,plsplot=False):
  #generate random single-spot lightcurve
  intensity, tparams, phase = genrandlc()
  print 'true params:', tparams
  print 'lightcurve:'
  print intensity
  print
  
  #fit using leastsq, starting at various points in param space
  opts,dists,p0s = lcfits(phase,intensity,tparams)
  
  #if no fits found, plot original curve and exit
  if not opts:
    print 'no fits found =('
    plt.plot(phase, intensity, 'r.')
    plt.xlabel('phase')
    plt.ylabel('intensity')
    plt.show()
    return opts,dists,p0s
  
  #find best fit
  bestdist, bdi = min((val, idx) for (idx, val) in enumerate(dists))
  bestparams = opts[bdi]
  bestp0s = p0s[bdi]
  bestlc = lcspot(phase,bestparams)
  print 'best fit:'
  print 'opt:',bestparams,'dist:',bestdist,'p0:',bestp0s
  print 'best lightcurve:'
  print bestlc
  print
  
  #print all fits
  if printall:
    print 'all fits:'
    for i in range(0,len(opts)):
      print 'opt:',opts[i],'dist:',dists[i],'p0:',p0s[i]
    print
  
  #plot
  if plsplot:
    plt.plot(phase, intensity, 'r.', phase, bestlc, 'b')
    plt.xlabel('phase')
    plt.ylabel('intensity')
    plt.show()
  
  return intensity, phase, tparams, opts, dists, p0s


def lcfit(phase,intensity):
  return leastsq(lcspotdiffs, [90.0,180.0,45.0,7.5], args=(phase,intensity))


def lcfits(phase,intensity,tps):
  opts = []
  dists = []
  p0s = []
  for inc in range(15,180,30):
   for lon in [180.0]:
    for lat in range(7,90,15):
     for rad in range(6,10,1):
      p0 = [inc,lon,lat,rad]
      try:
        fps, iem = leastsq(lcspotdiffs, p0, args=(phase,intensity))
        if iem in [1,2,3,4]:
          opts.append(fps)
          dists.append(sum(((fps/tps)-1)**2))
          p0s.append(p0)
      except RuntimeError as e:
        print "Runtime Error({0}): {1}".format(e.errno, e.strerror)
  return (opts,dists,p0s)


def main():
  if '-1' in sys.argv:
    randfit('-p' in sys.argv)
  else:
    randfit2(printall='-a' in sys.argv, plsplot='-p' in sys.argv)


main()
