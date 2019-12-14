

from matplotlib import pylab
from pylab import *

import numpy as np
import random


def column(matrix, i):
    return [row[i] for row in matrix]


tmax=1000 # Runtime in ms
N=1000 #no of neurons
Ne=800  #no of exc neurons             
Ni=200 #no of inh neurons      

w_upper = 2.0
w_lower = -2.0

#np.random.seed(1)

                    # Excitatory neurons    Inhibitory neurons
a = np.concatenate((np.full(Ne,0.02), np.full(Ni,0.01) ))
b = np.concatenate((np.full(Ne,0.2),  np.full(Ni,0.2) ))
c = np.concatenate((np.full(Ne,-65),  np.full(Ni,-65) ))
d = np.concatenate((np.full(Ne,8),    np.full(Ni,2) ))
#S = np.concatenate((0.5*np.random.uniform(0,1,(Ne,N)),    np.random.uniform(-1,0,(Ni,N)) )) # Excitatory neurons (positive weights only)    Inhibitory neurons (negative weights only) , multplying 0.5 for exh neurons seems to provide a more balanced NN 
Se = 0.5*np.random.rand(Ne+Ni,Ne)
Si = -np.random.rand(Ne+Ni,Ni)
S = np.hstack((Si,Se))
#S = S.T #flip the array so that each row has inh neruons at the end
np.savetxt("inital_weights.csv", S, delimiter=",")

print("initialized and saved inital weights")


#totzeros=0             # Count of zero-weight synapses
plasticity=19          # Plastic changes happen after the first 20 ms of runtime (in python 0-19)
v=-65*np.ones(N)    # Initial values of v
u=b * v                 # Initial values of u
firings=[]           # spike timings
LFP=np.zeros((tmax,3))
ADJ=np.zeros((tmax,N))
#count=0
#count2=0


#MAKE THE NETWORK SPARSELY CONNECTED 
M=np.random.rand(N,N) #Create an array of the given shape and populate it with random samples from a uniform distribution over [0, 1].
p=0.7                  # 30% of zero-weight synapses when p=0.7, this variable is used when the network is calculating spareness, anything greater than p is removed
for i in arange(tmax):
    for jjj in arange(tmax):
        if M[i,jjj]>p:
             S[i][jjj]=0
             #totzeros=totzeros+1
        
        if i==jjj:
             S[i][jjj]=0 #no self connections
             #totzeros=totzeros+1
SPARSE=S

np.savetxt("sparse_weights.csv", SPARSE, delimiter=",")

print("initialized and saved sparsed weights")

for t in arange(0,tmax): # simulation of runtime ms
    print("timestep: " + str(t))
    I=np.concatenate((5*np.random.randn(Ne),2*np.random.randn(Ni))) # thalamic input (apparently there is a random current firing between synapses, even during resting, this is seen in both reference examples), this is replaced by the input from spike trains and the rest should probably become zero, its only thalamic for randomized input
    fired=[i for i,x in enumerate(v) if x > 30]           # indices of spikes #spikes here are 30 not 35 (get fired neruons from previous step)

    for jj in arange(0,len(fired)):
        ADJ[t][fired[jj]]=1            # RECORD WHICH NEURONS FIRE AT EACH MS
 
    if t>plasticity:        #t>20 for STDP on, t>tmax for STDP off
        for j in arange(0,N):
            #count2=count2+1
            if ADJ[t][j]==1:   #neuron j fired at time t (for each neuron fired at that step)
                #count=count+1;
                for x in arange(0, plasticity):     #range of ms over which synaptic strengthening and weakening occur (go back 20 (plasticity) steps and look for fired neurons)
                    for k in arange(0,N): #iterate over all the neurons at that timestep looking for the neurons that fired
                        if ADJ[t-x][k]==1: # go backwards in timestep looking for the closest fired at 20 : 20,19,18,17

                            S[j][k]=S[j][k]*(1+(0.9*exp(x/20.0))) # synaptic weight of pair undergoes larger increment if dT is smaller and negative (increase weight i think?)
                            S[k][j]=S[k][j]*(1+(-0.925*exp(-x/20.0))) # synaptic weight of pair undergoes larger decrement is dT is smaller and positive (decrease weight i think?)

                            # set a maximum value for synaptic weights
                            if S[j][k]>w_upper:
                                S[j][k]=w_upper
                            
                            if S[j][k]<w_lower:
                                S[j][k]=w_lower
                            
                            if S[k][j]>w_upper:
                                S[k][j]=w_upper
                            
                            if S[k][j]<w_lower:
                                S[k][j]=w_lower

    #save the firings to print a chart
    print("fired: " + str(fired))
    firTemp = [[t,fired[i]] for i in arange(0,len(fired))]
    #print("firTemp: " + str(firTemp))
    if firTemp != []:
        firings = firings+ firTemp


    v[fired] = c[fired]
    u[fired]=u[fired]+d[fired]

    I=I+sum((S[:][fired]).T,1) #adjustment of curent based on synaptic weight value. get the columns of the fired neurons (meaning get all the post fired neuons, i think?), then sum the weights, and add to I (not 100% sure of the logic) 

    v=v+0.5*(0.04*v**2+5 * v+140-u+I) # step 0.5 ms
    v=v+0.5*(0.04*v**2+5 * v+140-u+I) # for numerical stability (apparently, fire twice per second)
    u=u+a*(b*v-u)                 # update u

    LFP[t][0]=sum(v[0:Ne])         # sum of voltages of excitatory per ms 
    LFP[t][1]=sum(v[Ne:N])     # sum of voltages of inhibitory per ms
    LFP[t][2] = sum(v)            # sum of all voltages for each ms

np.savetxt("final_weights.csv", S, delimiter=",")

figure()
tvec = arange(0, tmax)
fig,(vax1,vax2) = subplots(2)
vax1.plot(column(firings,0), column(firings,1), 'r.',label='Neurons',markersize=2)
xlabel('Time [ ms ]')
ylabel('Membrane voltage [mV]')
title( 'Spike raster - RS,FS STDP in a random network' )

vax2.plot(tvec, column(LFP,0), 'b', label = '# sum of voltages of excitatory per ms ')
vax2.plot(tvec, column(LFP,1), 'r', label = '# sum of voltages of inhibitory per ms')
vax2.plot(tvec, column(LFP,2), 'g', label = '# sum of all voltages for each ms')
xlabel('Time [ ms ]')
ylabel('Membrane voltage [mV]')
title( 'LFP - RS,FS STDP in a random network' )
show()




