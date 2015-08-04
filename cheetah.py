
from lcspotNEW import *
from lcsinglefitNEW import *
from lcmultifitNEW import *
from readfits import *

import time

def cheetah(fid):
  t,f,err = readfits('testdata/kplr'+str(fid)+'-2009350155506_llc.fits')
  lctuple = (t,f)
  t0 = time.time()
  paramsets = ratchetfit(lctuple, 3, p0s=[], initsteps=20, nclusters=30, threshratio=2, plsprint='some')
  t1 = time.time()
  print 'elapsed time:',(t1-t0)
  return paramsets
