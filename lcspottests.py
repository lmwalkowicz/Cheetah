from numpy import *
from lcspot import lcspot

jdlc1raw = loadtxt('jdlc1.txt')
jdlc1phase = jdlc1raw[:,0]
jdlc1expi = jdlc1raw[:,1]
jdlc1par = array([0.45, 0.3, 90.0, 180.0, 0.0, 4.0, 0.67])
jdlc1i = lcspot(jdlc1phase, jdlc1par)
print max(jdlc1i - jdlc1expi)

jdlc2raw = loadtxt('jdlc2.txt')
jdlc2phase = jdlc2raw[:,0]
jdlc2expi = jdlc2raw[:,1]
jdlc2par = array([0.45, 0.3, 90.0, 180.0, 20.0, 4.0, 0.67])
jdlc2i = lcspot(jdlc2phase, jdlc2par)
print max(jdlc2i - jdlc2expi)

jdlc4raw = loadtxt('jdlc4.txt')
jdlc4phase = jdlc4raw[:,0]
jdlc4expi = jdlc4raw[:,1]
jdlc4par = array([0.45, 0.3, 70.0, 180.0, 0.0, 4.0, 0.67])
jdlc4i = lcspot(jdlc4phase, jdlc4par)
print max(jdlc4i - jdlc4expi)
