

from matplotlib import pylab
from pylab import *

import numpy as np
import random

dt = 1
tmax=1000
N=1000 #no of neurons

a = 0.02
b = 0.2
c= -65
d = 8
S = 0
I = 10 #I is always applied continuous spike train current
I2 = 0
I3 = 0 #I3 only ever spikes if I1 or I2 spiked

v1 = np.zeros([tmax])
v2 = np.zeros([tmax])
v3 = np.zeros([tmax])

u1 = np.zeros([tmax]) 
u2 = np.zeros([tmax])
u3 = np.zeros([tmax])

dv1 = np.zeros([tmax])
dv2 = np.zeros([tmax])
dv3 = np.zeros([tmax])

fired_neurons = [[] for i in arange(0,1000)]

w_upper = 2.0
w_lower = -2.0

v1v3W = 2.0
v2v3W = 2.0

v1v3W_changes = np.full(tmax,2.0)
v2v3W_changes = np.full(tmax,2.0)

for t in arange(0,tmax-1):
    print("timestep: " + str(t))
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
        fired_neurons[t].append('v1')


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
        fired_neurons[t].append('v2')

    #check neuon 1 and 2 for spike at this point in time, if spiked, resolve action potiential of neuron 3
    #it is not a multiplication of weights like in traditional methods
    #for a given post neuron at that timestep, the weights of all synapses for fired pre neurons are added together
    #then added to the current 9again not multiplied)

    #first check i1 synapse
    added_weights=0
    if v1[t] == 35:
        added_weights = added_weights + v1v3W
        I3=10

    #first check i2 synapse and resolve    
    if v2[t] == 35:
        added_weights = added_weights + v2v3W
        I3=10

    #doing it the above way ensure I3 current only eveer spikes if v1 of v2 ever spiked
    #then add the weights to I3
    if I3==10:
        I3=I3+added_weights

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
        fired_neurons[t].append('v3')

    I3=0

    #add the minuimum 20ms check
    #do stdp
    if len(fired_neurons[t]) > 0:
        for fIndex in arange(0,len(fired_neurons[t])):

            #flipping this seems to create different results, check it later

            if(fired_neurons[t][fIndex]=='v1'):
                for lookBack in arange(1,20):
                    if 'v3' in fired_neurons[t-lookBack]:
                        v1v3W=v1v3W*(1+(-0.925*exp(-lookBack/20.0))) 
                        

            if(fired_neurons[t][fIndex]=='v2'):
                for lookBack in arange(1,20):
                    if 'v3' in fired_neurons[t-lookBack]:
                        v2v3W=v2v3W*(1+(-0.925*exp(-lookBack/20.0))) 

            if(fired_neurons[t][fIndex]=='v3'):
                for lookBack in arange(1,20):
                    if 'v1' in fired_neurons[t-lookBack]:
                        v1v3W=v1v3W*(1+(0.9*exp(lookBack/20.0))); 

                    if 'v2' in fired_neurons[t-lookBack]:
                        v2v3W=v2v3W*(1+(0.9*exp(lookBack/20.0))); 
                        
            if v1v3W > w_upper:
                v1v3W = w_upper

            if v1v3W < w_lower:
                v1v3W = w_lower

            if v2v3W > w_upper:
                v2v3W = w_upper

            if v2v3W < w_lower:
                v2v3W = w_lower

    v1v3W_changes[t] = v1v3W
    v2v3W_changes[t] = v2v3W




figure()
tvec = arange(0, tmax, dt)
fig,(vax1,vax2) = subplots(2)
vax1.plot(tvec, v1, 'r', label = 'v1 Voltage trace')
vax1.plot(tvec, v2, 'g', label = 'v2 Voltage trace')
vax1.plot(tvec, v3, 'b', label = 'v3 Voltage trace')
xlabel('Time [ ms ]')
ylabel('Membrane voltage [mV]')
title( 'voltage spikes' )

vax2.plot(tvec, v1v3W_changes, 'b', label = 'v1-v3 Weight Change')
vax2.plot(tvec, v2v3W_changes, 'g', label = 'v2-v3 Weight Change')
xlabel('Time [ ms ]')
ylabel('Weight')
title( 'weight changes' )
show()


    

