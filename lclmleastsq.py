
from numpy import *
from scipy.optimize import leastsq
from lcspot import *
import matplotlib.pyplot as plt
import sys

def randfit(plsplot=False):
  #generate random single-spot lightcurve
  phase = arange(0,1,.025)
  tinc = random.uniform(0,180)
  tlon = random.uniform(0,360)
  tlat = random.uniform(0,90)
  trad = random.uniform(5,10)
  tparams = array([tinc, tlon, tlat, trad])
  i = lcspot(phase, tparams)
  print 'true params:', tparams
  print 'lightcurve:', i
  
  #fit using leastsq, starting in the center of param space 
  fitparams, iem = lcfit(phase,i)
  print 'fit params:', fitparams
  print 'iem:', iem
  
  #plot
  if plot:
    if iem in [1,2,3,4]:
      plt.plot(phase, i, 'r', phase, lcspot(phase,fitparams), 'b')
      plt.xlabel('phase')
      plt.ylabel('intensity')
      plt.show()
    else:
      plt.plot(phase, i, 'r', phase)
      plt.xlabel('phase')
      plt.ylabel('intensity')
      plt.show()
  return fitparams, iem


def randfit2(printall=False,plsplot=False):
  #generate random single-spot lightcurve
  phase = arange(0,1,.025)
  tinc = random.uniform(0,180)
  tlon = random.uniform(0,360)
  tlat = random.uniform(0,90)
  trad = random.uniform(5,10)
  tparams = array([tinc, tlon, tlat, trad])
  i = lcspot(phase, tparams)
  print 'true params:', tparams
  print 'lightcurve:'
  print i
  print
  
  #fit using leastsq, starting at various points in param space
  opts,dists,p0s = lcfits(phase,i,tparams)
  
  #if no fits found, plot original curve and exit
  if not opts:
    print 'no fits found =('
    plt.plot(phase, i, 'r')
    plt.xlabel('phase')
    plt.ylabel('intensity')
    plt.show()
    return opts,dists,p0s
  
  #find best fit
  mindist, mdi = min((val, idx) for (idx, val) in enumerate(dists))
  print 'best fit:'
  print 'opt:',opts[mdi],'dist:',dists[mdi],'p0:',p0s[mdi]
  print 'best lightcurve:'
  print lcspot(phase,opts[mdi])
  print
  
  #print all fits
  if printall:
    print 'all fits:'
    for i in range(0,len(opts)):
      print 'opt:',opts[i],'dist:',dists[i],'p0:',p0s[i]
    print
  
  #plot
  if plsplot:
    plt.plot(phase, i, 'r', phase, lcspot(phase,opts[mdi]), 'b')
    plt.xlabel('phase')
    plt.ylabel('intensity')
    plt.show()
  
  return opts,dists,p0s


def lcfit(phase,i):
  return leastsq(lcspotdiffs, [90.0,180.0,45.0,7.5], args=(phase,i))


def lcfits(phase,i,tps):
  opts = []
  dists = []
  p0s = []
  for inc in range(15,180,30):
   for lon in [180.0]:
    for lat in range(7,90,15):
     for rad in range(6,10,1):
      p0 = [inc,lon,lat,rad]
      try:
        fps, iem = leastsq(lcspotdiffs, p0, args=(phase,i))
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
