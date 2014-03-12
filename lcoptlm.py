
from numpy import *
from scipy.optimize import curve_fit
from lcspot2 import lcspot2

def randfit():
  phase = arange(0,1,.025)
  tinc = random.uniform(0,180)
  #tlon = random.uniform(0,360)
  tlon = 180.0
  tlat = random.uniform(0,90)
  trad = random.uniform(5,10)
  i = lcspot2(phase, tinc, tlon, tlat, trad)
  tparams = array([tinc, tlon, tlat, trad])
  print 'true params:', tparams
  print 'lightcurve:', i
  return lcfits(phase,i,tparams)

def lcfits(phase,i,tparams):
  opts = []
  dists = []
  p0s = []
  for inc in range(15,180,30):
   for lon in [180.0]:
    for lat in range(7,90,15):
     for rad in range(6,10,1):
      p0 = [inc,lon,lat,rad]
      try:
        popt, pcov = curve_fit(lcspot2, phase, i, p0)
        opts.append(popt)
        dists.append(sum(((popt/tparams)-1)**2))
        p0s.append(p0)
      except RuntimeError:
        print 'could not find opt for p0 =', p0
  return (opts,dists,p0s)

