
# port of lcspot.pro
# Analytic model of intensity variations, due to a circular spot on a
# rotating stars, based on equations in Eker (1994, ApJ, 420, 373).
#
# Inputs:
#   phase (array[nphi]) rotational phases for output light curve, where
#     phase=0.0 is when stellar longitude 0 crosses disk center.
#   par (array[7]) stellar and spot parameters
#     par[0:6] = [limb1, limb2, inc_deg, lon_deg, lat_deg, rad_deg, iratio]
#     limb1: linear term in limb-darkening law
#       I[mu] = I[0] .* (1 - limb1.*(1-mu) - limb2.*(1-mu.^2)), where mu=cos(theta)
#         Form used by Claret (2004)#   Eker (1994) uses u1=limb1, u2=-limb2.
#         Limb-darkening coefficients depend on filter bandpass.
#     limb2: quadratic term in limb-darkening law [see limb1]
#     inc_deg: inclination of stellar rotation axis (in degrees)
#     lon_deg: longitude of spot center (in degrees)
#     lat_deg: latitude of spot center (in degrees)
#     rad_deg: angle radius of spot (in degrees)
#     iratio: spot intensity divided by non-spot intensity
#       0:black, <1:dark, 1:photosphere, >1:bright
#

from numpy import *
import itertools

def lcspot(time, params):
  #Internal constants
  d2r = pi / 180
  
  #Impose Boundaries
  bps = boundparams(params)
  
  #Extract stellar and spot properties from parameter vector
  inc_deg = bps[0]
  lon_deg = bps[1]
  lat_deg = bps[2]
  rad_deg = bps[3]
  period = bps[4]
  limb1 = 0.45    #linear coefficient in limb-darkening law
  limb2 = 0.3     #quadratic coefficient in limb-darkening law
  iratio = 0.67   #intensity in spot / intensity out of spot
  
  #Convert input angles from degrees to radians
  inc = inc_deg * d2r                  #stellar inclination
  lam = lon_deg * d2r                  #spot longitude
  bet = lat_deg * d2r                  #spot latitude
  rad = rad_deg * d2r                  #spot radius
  
  #Useful scalar quantities
  cosrad = cos(rad)
  sinrad = sin(rad)
  
  #Calculate a vector of rotational phases (could be passed as an argument)
  phase = time / period
  phi = 2.0 * pi * phase
  nphi = len(phi)
  
  #Calculate angle "theta" between two vector originating at spot center:
  # Vector 1) normal to stellar surface directed away from center of star
  # Vector 2) directed towards the observer.
  #See equation 20 of Eker (1994).
  costhe0 = cos(inc) * sin(bet) + sin(inc) * cos(bet) * cos(phi-lam)
  
  #Useful quantities
  sinthe0 = sqrt(1.0 - costhe0**2)
  the0 = arccos(costhe0)
  
  #Find rotational phases when spot is full, gibbous, crescent, occulted
  jf = flatnonzero(the0 <= pi/2-rad)
  nf = len(jf)
  jg = flatnonzero(logical_and(the0 > pi/2-rad, the0 <= pi/2))
  ng = len(jg)
  jc = flatnonzero(logical_and(the0 > pi/2, the0 <= pi/2+rad))
  nc = len(jc)
  jo = flatnonzero(the0 > pi/2+rad)
  no = len(jo)
  
  #Allocate vectors for intensity integrals.
  ic = zeros(nphi)    #constant intensity term
  il = zeros(nphi)    #linear intensity term
  iq = zeros(nphi)    #quadratic intensity term
  
  #
  # Rotational phases when spot is full (entirely visible)
  #
  
  #Useful quantities for full phases.
  if nf >= 1:
    #unused: the0_f = the0[jf]        #angle between spot and LOS
    costhe0_f = costhe0[jf]
    sinthe0_f = sinthe0[jf]
    
    #Calculate intensity integrals for phases when spot is fully visible
    #See equations 13a, 13b, and 13c of Eker (1994)
    ic[jf] = pi * sin(rad)**2 * costhe0_f
    il[jf] = 2*pi/3 * (1 - cosrad**3) - pi * cosrad * sinrad**2 * sinthe0_f**2
    iq[jf] = pi/2 * (1 - cosrad**4) * costhe0_f**3 + 3*pi/4 * sinrad**4 * costhe0_f * sinthe0_f**2
    
  #
  # Rotational phases when spot is gibbous (more than half visible).
  #
  
  #Useful quantities for gibbous phases.
  if ng >= 1:
    the0_g = the0[jg]                #angle between spot and LOS
    costhe0_g = costhe0[jg]
    sinthe0_g = sinthe0[jg]
    
    #Calculate integration limits for partially visible spot.
    #See equation 16 of Eker (1994).
    cosphi0_g = - 1.0 / ( tan(the0_g) * tan(rad) )
    rad0_g = abs( the0_g - pi/2 )
    
    #Useful quantities.
    #if any(abs(cosphi0_g) > 1.0): print 'IVW @ 114 & 115',params
    phi0_g = arccos(cosphi0_g)
    sinphi0_g = sqrt(1.0 - cosphi0_g**2)
    cosrad0_g = cos(rad0_g)
    sinrad0_g = sin(rad0_g)
    
    #Auxiliary quantities for gibbous phases.
    #See unnumbered equations that follow equation 18b in Eker (1994).
    k1_g = ((pi - phi0_g) / 4) * (cosrad0_g**4 - cosrad**4)
    k2_g = (sinphi0_g / 8) * ( rad0_g - rad + 0.5 * ( sin(2*rad)  * cos(2*rad) - sin(2*rad0_g) * cos(2*rad0_g) ) )
    k3_g = (1.0 / 8) * (pi - phi0_g - sinphi0_g * cosphi0_g) * (sinrad**4 - sinrad0_g**4)
    k4_g = - (sinphi0_g - sinphi0_g**3 / 3) * ( (3.0 / 8) * (rad - rad0_g) + (1.0 / 16) * ( sin(2*rad)  * (cos(2*rad)  - 4) - sin(2*rad0_g) * (cos(2*rad0_g) - 4) ) )
    
    #Corrections to intensity integrals for gibbous phases.
    #See equations 18a and 18b in Eker (1994).
    cl_g = ((pi - phi0_g) / 3) * (cosrad**3 - cosrad0_g**3) * (1 - 3*costhe0_g**2) - (pi - phi0_g - sinphi0_g * cosphi0_g) * (cosrad - cosrad0_g) * sinthe0_g**2 - (4.0 / 3) * sinphi0_g * (sinrad**3 - sinrad0_g**3) * sinthe0_g * costhe0_g - (1.0 / 3) * sinphi0_g * cosphi0_g * (cosrad**3 - cosrad0_g**3) * sinthe0_g**2
    cq_g = 2 * costhe0_g**3 * k1_g + 6 * costhe0_g**2 * sinthe0_g * k2_g + 6 * costhe0_g * sinthe0_g**2 * k3_g + 2 * sinthe0_g**3 * k4_g
    
    #Calculate intensity integrals for phases when spot is fully visible.
    #Constant intensity integral. See equation 17 of Eker (1994)
    #if any(abs(cosrad / sinthe0_g) > 1.0): print 'IVW @ 134',params
    ic[jg] = phi0_g * costhe0_g * sinrad**2 - arcsin(cosrad / sinthe0_g) - 0.5 * sinthe0_g * sinphi0_g * sin(2*rad) + pi/2
    
    #Apply corrections to linear and quadratic intensity integrals.
    il[jg] = 2*pi/3 * (1 - cosrad**3) - pi * cosrad * sinrad**2 * sinthe0_g**2 - cl_g
    iq[jg] = pi/2 * (1 - cosrad**4) * costhe0_g**3 + 3*pi/4 * sinrad**4 * costhe0_g * sinthe0_g**2 - cq_g
    
  #
  # Rotational phases when spot is crescent (less than half visible).
  #
  
  #Useful quantities for crescent phases.
  if nc >= 1:
    the0_c = the0[jc]        #angle between spot and LOS
    costhe0_c = costhe0[jc]
    sinthe0_c = sinthe0[jc]
    
    #Calculate integration limits for partially visible spot.
    #See equation 16 of Eker (1994).
    cosphi0_c = - 1.0 / ( tan(the0_c) * tan(rad) )
    rad0_c = abs( the0_c - pi/2 )
    
    #Useful quantities.
    #if any(abs(cosphi0_c) > 1.0): print 'IVW @ 157 & 158',params
    phi0_c = arccos(cosphi0_c)
    sinphi0_c = sqrt(1.0 - cosphi0_c**2)
    cosrad0_c = cos(rad0_c)
    sinrad0_c = sin(rad0_c)
    
    #Auxiliary quantities for crescent phases.
    #See unnumbered equations that follow equation 18b in Eker (1994).
    k1_c = (phi0_c / 4) * (cosrad0_c**4 - cosrad**4)
    k2_c = - (sinphi0_c / 8) * ( rad0_c - rad + 0.5 * ( sin(2*rad)  * cos(2*rad) - sin(2*rad0_c) * cos(2*rad0_c) ) )
    k3_c = (1.0 / 8) * (phi0_c + sinphi0_c * cosphi0_c) * (sinrad**4 - sinrad0_c**4)
    k4_c = (sinphi0_c - sinphi0_c**3 / 3) * ( (3.0 / 8) * (rad - rad0_c) + (1.0 / 16) * ( sin(2*rad)  * (cos(2*rad)  - 4) - sin(2*rad0_c) * (cos(2*rad0_c) - 4) ) )
    
    #Corrections to intensity integrals for crescent phases.
    #See equations 18a and 18b in Eker (1994).
    #unused: cl_c = ((pi - phi0_c) / 3) * (cosrad**3 - cosrad0_c**3) * (1 - 3*costhe0_c**2) - (pi - phi0_c - sinphi0_c * cosphi0_c) * (cosrad - cosrad0_c) * sinthe0_c**2 - (4.0 / 3) * sinphi0_c * (sinrad**3 - sinrad0_c**3) * sinthe0_c * costhe0_c - (1.0 / 3) * sinphi0_c * cosphi0_c * (cosrad**3 - cosrad0_c**3) * sinthe0**2
    cq_c = 2 * costhe0_c**3 * k1_c + 6 * costhe0_c**2 * sinthe0_c * k2_c + 6 * costhe0_c * sinthe0_c**2 * k3_c + 2 * sinthe0_c**3 * k4_c
    
    #Calculate intensity integrals for phases when spot is fully visible.
    #Constant intensity integral. See equation 17 of Eker (1994)
    #if any(abs(cosrad / sinthe0_c) > 1.0): print 'IVW @ 177',params
    ic[jc] = phi0_c * costhe0_c * sinrad**2 - arcsin(cosrad / sinthe0_c) - 0.5 * sinthe0_c * sinphi0_c * sin(2*rad) + pi/2
    
    #Linear intensity integral. See equation 19a of Eker (1994)
    il[jc] = (phi0_c / 3) * (cosrad**3 - cosrad0_c**3) * (1 - 3 * costhe0_c**2) - (phi0_c + sinphi0_c * cosphi0_c) * (cosrad - cosrad0_c) * sinthe0_c**2 + (4.0 / 3) * sinphi0_c * (sinrad**3 - sinrad0_c**3) * sinthe0_c * costhe0_c + (1.0 / 3) * sinphi0_c * cosphi0_c * (cosrad**3 - cosrad0_c**3) * sinthe0_c**2
    
    #Apply corrections to linear and quadratic intensity integrals.
    iq[jc] = cq_c
    
  #
  # Rotational phases when spot is completely occulted (back of star).
  #
  
  if no >=1:
    ic[jo] = 0.0
    il[jo] = 0.0
    iq[jo] = 0.0
  
  #
  # Calculate light curve. Equation 12c from Eker (1994).
  #
  
  lc = 1.0 + (iratio - 1.0) / (pi * (1.0 - limb1/3.0 + limb2/6.0)) * ((1.0 - limb1 + limb2)*ic + (limb1 - 2.0 * limb2)*il + limb2*iq)
  
  return lc


### Single Spot Helper Functions ###

def lcspotdiffs(params, time, data):
  return data - lcspot(time, params)

def lcspotdiffs2(params, time, data):
  diffs = data - lcspot(time, params)
  penalty = 0.0
  if params[3] < 0.01:
    penalty = penalty + (((-params[3] + 0.01)/20.0)**2)/len(diffs)
  if params[4] < 0.01:
    penalty = penalty + (((-params[4] + 0.01)/20.0)**2)/len(diffs)
  return sqrt(diffs**2 + penalty)

def lcspotdiffs3(params, time, data):
  diffs = data - lcspot(time, params)
  posdiffs = diffs[diffs > 0.0]
  negdiffs = diffs[diffs <= 0.0]
  posdiffs = 4*posdiffs
  diffs2 = concatenate((posdiffs, negdiffs))
  penalty = 0.0
  if params[3] < 0.01:
    penalty = penalty + (((-params[3] + 0.01)/20.0)**2)/len(diffs)
  if params[4] < 0.01:
    penalty = penalty + (((-params[4] + 0.01)/20.0)**2)/len(diffs)
  return sqrt(diffs2**2 + penalty)

def lcspotsse(params, time, data):
  return sum(lcspotdiffs(params, time, data)**2)

def lcspotsse2(params, time, data):
  return sum(lcspotdiffs2(params, time, data)**2)

def lcspotsse3(params, time, data):
  return sum(lcspotdiffs3(params, time, data)**2)


### Fixed Inclination Helper Functions ###

def lcspotfi(time, inc, sparams):
  return lcspot(time, [inc] + list(sparams))

def lcspotdiffsfi2(params, time, data, inc):
  return lcspotdiffs2(array([inc] + list(params)), time, data)

def lcspotdiffsfi3(params, time, data, inc):
  return lcspotdiffs3(array([inc] + list(params)), time, data)

def lcspotssefi(params, time, data, inc):
  return lcspotsse(array([inc] + list(params)), time, data)

def lcspotssefi2(params, time, data, inc):
  return sum(lcspotdiffsfi2(params, time, data, inc)**2)

def lcspotssefi3(params, time, data, inc):
  return sum(lcspotdiffsfi3(params, time, data, inc)**2)


### Fixed Star Helper Functions ###

def lcspotfstar(time, inc, teq, alpha, sparams):
  lat = sparams[1]
  period = teq / (1 - alpha * sin(lat)**2)
  return lcspot(time, [inc] + list(sparams) + [period])

def lcspotdiffsfstar2(params, time, data, inc, teq, alpha):
  lat = params[1]
  period = teq / (1 - alpha * sin(lat)**2)
  return lcspotdiffs2(array([inc] + list(params) + [period]), time, data)

def lcspotdiffsfstar3(params, time, data, inc, teq, alpha):
  lat = params[1]
  period = teq / (1 - alpha * sin(lat)**2)
  return lcspotdiffs3(array([inc] + list(params) + [period]), time, data)

def lcspotssefstar(params, time, data, inc, teq, alpha):
  return sum(data - lcspotfstar(time, inc, teq, alpha, params))**2


### Multi-Spot Helper Functions ###

def lccomb(lc1, lc2):
  return lc1 + lc2 - 1

def lcsep(combined, component):
  return combined - component + 1

def lcmultispot(time, pset):
  pset = boundparamsms(pset)
  inc = pset[0]
  teq = pset[1]
  alpha = pset[2]
  spots = pset[3:]
  intensity = 1
  for spot in spots:
    tspot = teq / (1.0 - alpha * sin(spot[1] * pi / 180.0)**2)
    intensity = lccomb(intensity, lcspot(time,[inc]+list(spot)+[tspot]))
  return intensity

def lcmultispotdiffs(params,time,data):
  gpl = fpl2gpl(params)
  diffs = abs(data - lcmultispot(time, gpl))
  penalty = 0.0
  for s in gpl[3:]:
    if s[2] < 0.01:
      penalty = penalty + (((-s[2] + 0.01)/20.0)**2)/len(diffs)
  return sqrt(diffs**2 + penalty)

def lcmultispotsse(params, time, data):
  diffs = abs(data - lcmultispot(time, fpl2gpl(params)))
  return sum(diffs**2)


### Other Helper Functions ###

def paramdist(fps, tps, scalevals=array([45.0, 180.0, 90.0, 17.5, 25.0])):
  return sqrt(sum(((array(fps) - array(tps))/array(scalevals))**2))

def spotparamdist(fps, tps, scalevals=array([180.0, 90.0, 17.5])):
  return paramdist(fps, tps, scalevals)

def multispotparamdist(fps, tps, scalevals=array([180.0, 90.0, 17.5])):
  return sqrt(sum(array([paramdist(fps[i], tps[i], scalevals) for i in range(len(fps))])**2))
  
def spotmatch(fspots,tspots,scalevals=array([180.0, 90.0, 17.5])):
  nspots = len(fspots)
  if nspots == 1:
    return fspots, tspots, spotparamdist(fspots[0],tspots[0],scalevals)
  else:
    mindist = double('inf')
    for r in range(nspots):
      rfspots, rtspots, rdist = spotmatch(fspots[1:], tspots[:r]+tspots[r+1:], scalevals)
      totaldist = spotparamdist(fspots[0],tspots[r],scalevals)
      if totaldist < mindist:
        mindist = totaldist
        nfspots = fspots[0:1] + rfspots
        ntspots = tspots[r:r+1] + rtspots
    return nfspots, ntspots, mindist

def fpl2gpl(fpl):
  return [fpl[0], fpl[1], fpl[2]] + [[fpl[3+3*i],fpl[3+3*i+1],fpl[3+3*i+2]] for i in range((len(fpl)-3)/3)]

def gpl2fpl(gpl):
  return [gpl[0], gpl[1], gpl[2]] + [sp for sps in gpl[3:] for sp in sps]

def spacedvals(min, max, nvals):
  spacing = (max - min)/double(nvals)
  return arange(min+(spacing/2.0),max,spacing)

def spacedparams(nvals, minr=10, maxr=25, mint=1.0, maxt=45.0):
  incvals = spacedvals(0,90,nvals)
  lonvals = spacedvals(0,360,nvals)
  latvals = spacedvals(-90,90,nvals)
  radvals = spacedvals(minr,maxr,nvals)
  pervals = spacedvals(mint,maxt,nvals)
  return [list(x) for x in itertools.product(incvals,lonvals,latvals,radvals,pervals)]

def spacedparamsfi(nvals, minr=10, maxr=25, mint=1.0, maxt=45.0):
  lonvals = spacedvals(0,360,nvals)
  latvals = spacedvals(-90,90,nvals)
  radvals = spacedvals(minr,maxr,nvals)
  pervals = spacedvals(mint,maxt,nvals)
  return [list(x) for x in itertools.product(lonvals,latvals,radvals,pervals)]

def spacedparamsspot(nvals, minr=10, maxr=25):
  lonvals = spacedvals(0,360,nvals)
  latvals = spacedvals(-90,90,nvals)
  radvals = spacedvals(minr,maxr,nvals)
  return [list(x) for x in itertools.product(lonvals,latvals,radvals)]

def canonicize(inc,lon,lat):
  clat = ((lat + 180.0) % 360.0) - 180.0
  clon = lon
  if clat > 90.0:
    clat = 180.0 - clat
    clon = clon + 180.0
  if clat < -90.0:
    clat = -180.0 - clat
    clon = clon + 180.0
  cinc = ((inc + 180.0) % 360.0) - 180.0
  if cinc < 0.0:
    cinc = -cinc
    clon = clon + 180.0
  if cinc > 90.0:
    cinc = 180.0 - cinc
    clat = -clat
  clon = clon % 360.0
  return cinc,clon,clat

def boundparams(params):
  inc,lon,lat = canonicize(params[0],params[1],params[2])
  rad = (params[3] if params[3] >= 0.01 else 0.01)
  period = (params[4] if params[4] >= 0.01 else 0.01)
  return array([inc,lon,lat,rad,period])

def boundparamsfi(params):
  inc,lon,lat = canonicize(45.0,params[0],params[1])
  rad = (params[2] if params[2] >= 0.01 else 0.01)
  period = (params[3] if params[3] >= 0.01 else 0.01)
  return array([lon,lat,rad,period])

def boundparamsfstar(params):
  inc,lon,lat = canonicize(45.0,params[0],params[1])
  rad = (params[2] if params[2] >= 0.01 else 0.01)
  return array([lon,lat,rad])

def boundparamsms(gpl):
  inc = gpl[0]
  cinc = ((inc + 180.0) % 360.0) - 180.0
  if cinc < 0.0: cinc = -cinc
  if cinc > 90.0: cinc = 180.0 - cinc
  return [cinc,gpl[1],gpl[2]] + [list(canonicize(inc,s[0],s[1])[1:])+[(s[2] if s[2] >= 0.01 else 0.01)] for s in gpl[3:]]



