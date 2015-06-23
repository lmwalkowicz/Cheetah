
from numpy import *
import matplotlib.pyplot as plt
from scipy.interpolate import *
plt.ion()

def splinefit(x,y,k=3,s=.005):
  uspl = UnivariateSpline(x,y,k=k,s=s)
  r = arange(x[0],x[-1],.01)
  n = uspl.get_knots()
  plt.figure()
  plt.plot(x,y,'r.',r,uspl(r),'b',n,uspl(n),'o')
  return uspl

def lsqsplinefit(x,y,k=3,n=100):
  t = arange(x[0],x[-1],(x[-1]-x[0])/(n+1))[1:]
  spl = LSQUnivariateSpline(x,y,t,k=k)
  r = arange(x[0],x[-1],.01)
  n = spl.get_knots()
  plt.figure()
  plt.plot(x,y,'r.',r,spl(r),'b',n,spl(n),'o')
  return spl

def findoutliers(x,y):
  uspl1 = splinefit(x,y,s=.5)
  return uspl