
from numpy import isfinite
import pyfits

def readfits(file):
  lc = pyfits.getdata(file)
  t = lc.field('TIME')
  f = lc.field('PDCSAP_FLUX')
  err = lc.field('PDCSAP_FLUX_ERR')
  f = f[isfinite(t)]
  err = err[isfinite(t)]
  t = t[isfinite(t)]
  t = t[isfinite(f)]
  err = err[isfinite(f)]
  f = f[isfinite(f)]
  f = f[isfinite(err)]
  t = t[isfinite(err)]
  err = err[isfinite(err)]
  sf = sorted(f)
  pseudomax = sf[int(len(sf)*.98)]
  f = f/pseudomax
  t = t - t[0]
  return t,f,err

