
from numpy import *
from scipy.optimize import curve_fit
from lcspot3 import lcspot3

def randfit():
  phase = arange(0,1,.025)
  #tinc = random.uniform(0,180)
  #tlon = random.uniform(0,360)
  #tlon = 180.0
  tlat = random.uniform(0,90)
  #trad = random.uniform(5,10)
  i = lcspot3(phase, tlat)
  tparams = array([tlat])
  print 'true params:', tparams
  #print 'lightcurve:', i
  return lcfits(phase,i,tparams)

def lcfits(phase,i,tparams):
  opts = []
  dists = []
  p0s = []
  for lat in range(7,90,15):
      p0 = [lat]
      #try:
      popt, pcov = curve_fit(lcspot3, phase, i, p0)
      opts.append(popt)
      dists.append(sum(((popt/tparams)-1)**2))
      p0s.append(p0)
      #except RuntimeError:
      #  print 'could not find opt for p0 =', p0
  return (opts,dists,p0s)

