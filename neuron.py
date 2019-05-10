#Exercise 1
#Simulate one excitatory neuron for 1000ms and plot the resulting voltage
#trace. Apply a current step (Iapp = 7 (picoAmperes)) between time 200ms and 700ms.
#https://en.wikipedia.org/wiki/Current_density  Current densities are usually expressed in pA⋅pF−1
from matplotlib import pylab
from pylab import *

# 1 ) initialize parameters
tmax = 1000
dt = 0.5

# 1 . 1 ) Neuron / Network parameters (related to last slide of introduction)
a = 0.02 #RS , IB: 0.02 , FS: 0.1
b = 0.2 #RS , IB , FS: 0.2
c = -65 #RS , FS: −65 IB: −55
d = 8 #RS: 8 , IB: 4 , FS: 2

# 1 . 2 ) Input Parameters
Iapp=10 #(this should be 7?)
tr= array([200,700] )/dt # stm time

# 2 ) reserve memory
T = int(ceil(tmax/dt))
v = zeros(T)
u = zeros(T)
v[0] = -70 # resting potential
u[0] = -14 # steady state

# 3 ) for−loop over time
for t in arange(T-1):
# 3 . 1 ) get input
    if t>tr[0] and t<tr[1] :
        I = Iapp
    else :
        I = 0

    if v[t] < 35:
# 3 . 2 ) update ODE (ordinary differential equations)
        dv = (0.04 * v[t] + 5) * v[t] + 140 - u[t]
        v[t + 1] = v[t] + (dv + I) * dt
        du = a * (b * v [t] - u[t])
        u[t + 1] = u[t] + dt * du
    else :
# 3 . 3 ) spike!
        v[t] = 35
        v[t + 1] = c
        u[t + 1] = u [t] + d


# 4 ) plot voltage trace
figure()
tvec = arange(0, tmax, dt)
plot(tvec, v, 'b', label = 'Voltage trace')
xlabel('Time [ ms ]')
ylabel('Membrane voltage [mV]')
title( 'A single qIF neuron with current step input6' )
show()