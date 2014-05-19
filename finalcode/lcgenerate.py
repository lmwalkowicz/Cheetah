
from numpy import *
from lcspot import *

def genrandlc(phase=arange(0,1,.025), noisefactor=0.02):
  while True:
    inc = random.uniform(0,90)
    lon = random.uniform(0,360)
    lat = random.uniform(-90,90)
    rad = random.uniform(10,25)
    params = array([inc, lon, lat, rad])
    intensity = lcspot(phase, params)
    if min(intensity) < .999: break
  if noisefactor > 0.0:
    lcrange = max(intensity) - min(intensity)
    noise = random.normal(0,lcrange*noisefactor,len(intensity))
    intensity = intensity + noise
  return phase, intensity, params


def genfinclc(tinc, phase=arange(0,1,.025), noisefactor=0.02):
  while True:
    lon = random.uniform(0,360)
    lat = random.uniform(-90,90)
    rad = random.uniform(10,25)
    params = array([inc, lon, lat, rad])
    intensity = lcspot(phase, params)
    if min(intensity) < .999: break
  if noisefactor > 0.0:
    lcrange = max(intensity) - min(intensity)
    noise = random.normal(0,lcrange*noisefactor,len(intensity))
    intensity = intensity + noise
  return phase, intensity, params


def genmultilc(phase=arange(0,1,.0125), noisefactor=0.02, nspots=2):
  while True:
    inc = random.uniform(0,90)
    lon = random.uniform(0,360)
    lat = random.uniform(-90,90)
    rad = random.uniform(10,25)
    intensity = lcspot(phase, [inc, lon, lat, rad])
    if min(intensity) < .999: break
  spots = [[lon, lat, rad]]
  pang = 90 - lat
  d2r = pi / 180
  spotsc = [[sin(d2r*pang)*cos(d2r*lon), sin(d2r*pang)*sin(d2r*lon), cos(d2r*pang), rad]]
  
  for i in range(1,nspots):
    while True:
      lon = random.uniform(0,360)
      lat = random.uniform(-90,90)
      rad = random.uniform(10,25)
      pang = 90 - lat
      x = sin(d2r*pang)*cos(d2r*lon)
      y = sin(d2r*pang)*sin(d2r*lon)
      z = cos(d2r*pang)
      collision = False
      for spotc in spotsc:
        angle = arccos(x*spotc[0]+y*spotc[1]+z*spotc[2])/d2r
        if angle < (rad + spotc[3]):
          collision = True
          break
      if not collision:
        intensity2 = lcspot(phase, [inc, lon, lat, rad])
        if min(intensity2) < .999: break
    intensity = lccomb(intensity, intensity2)
    spots.append([lon, lat, rad])
    spotsc.append([x, y, z, rad])
  if noisefactor > 0.0:
    lcrange = max(intensity) - min(intensity)
    noise = random.normal(0,lcrange*noisefactor,len(intensity))
    intensity = intensity + noise
  return phase, intensity, [inc] + spots

