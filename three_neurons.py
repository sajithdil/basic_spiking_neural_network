# 3 neurons
# 1st continuous spike
# 2nd random spike
# 3rd action potenital increases if 1 or 2 fires 

from matplotlib import pylab
from pylab import *

import numpy as np
import random

dt = 1
tmax=1000

a = 0.02
b = 0.2
c= -65
d = 8
S = 0
I = 10 #I is always applied continuous spike train current
I2 = 0
I3 = 0 #I3 only ever spikes if I1 or I2 spiked

v1 = np.zeros([1000])
v2 = np.zeros([1000])
v3 = np.zeros([1000])

u1 = np.zeros([1000]) 
u2 = np.zeros([1000])
u3 = np.zeros([1000])

dv1 = np.zeros([1000])
dv2 = np.zeros([1000])
dv3 = np.zeros([1000])

fired_neurons = []

for t in arange(0,tmax-1):
    s = random.randint(0,1)
    if s == 0:
        I2=0
    else:
        I2=10

    if v1[t]<35:
        dv = (0.04 * v1[t] + 5) * v1[t] + 140 - u1[t]
        dv1[t]=dv
        v1[t + 1] = v1[t] + (dv + I) * dt
        du = a * (b * v1[t] - u1[t])
        u1[t + 1] = u1[t] + dt * du 
    else:
        v1[t] = 35
        v1[t + 1] = c
        u1[t + 1] = u1[t] + d


    if v2[t]<35:
        dv = (0.04 * v2[t] + 5) * v2[t] + 140 - u2[t]
        dv2[t]=dv
        v2[t + 1] = v2[t] + (dv + I2) * dt
        du = a * (b * v2[t] - u2[t])
        u2[t + 1] = u2[t] + dt * du 
    else:
        v2[t] = 35
        v2[t + 1] = c
        u2[t + 1] = u2[t] + d

    #check neuon 1 and 2 for spike at this point in time, if spiked, resolve action potiential of neuron 3

    #first check i1 synapse and resolve
    if v1[t] == 35:
        I3=10
    else:
        I3=0

    if v3[t]<35:
        dv = (0.04 * v3[t] + 5) * v3[t] + 140 - u3[t]
        dv3[t]=dv
        v3[t + 1] = v3[t] + (dv + I3) * dt
        du = a * (b * v3[t] - u3[t])
        u3[t + 1] = u3[t] + dt * du 
    else:
        v3[t] = 35
        v3[t + 1] = c
        u3[t + 1] = u3[t] + d

    #first check i2 synapse and resolve    
    if v2[t] == 35:
        I3=10
    else:
        I3=0

    if v3[t]<35:
        dv = (0.04 * v3[t] + 5) * v3[t] + 140 - u3[t]
        dv3[t]=dv
        v3[t + 1] = v3[t] + (dv + I3) * dt
        du = a * (b * v3[t] - u3[t])
        u3[t + 1] = u3[t] + dt * du 
    else:
        v3[t] = 35
        v3[t + 1] = c
        u3[t + 1] = u3[t] + d


figure()
tvec = arange(0, tmax, dt)
plot(tvec, v1, 'r', label = 'Voltage trace')
plot(tvec, v2, 'g', label = 'Voltage trace')
plot(tvec, v3, 'b', label = 'Voltage trace')
xlabel('Time [ ms ]')
ylabel('Membrane voltage [mV]')
title( 'v3 voltage spikes' )
show()

    

