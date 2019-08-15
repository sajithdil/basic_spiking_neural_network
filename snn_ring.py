#Exercise 3
#Simulate 1000 neurons for 1000 ms and plot the resulting spikes. Assume
#that each neuron receives (random) 10% of the 100 Poisson spike trains of
#rate frate = 2 Hz between time 200 ms and 700 ms. Note that the neurons
#are not yet inter-connected.

from matplotlib import pylab
from pylab import *
from scipy.sparse import csr_matrix
from scipy.linalg import circulant

# 1 ) initialize parameters
tmax = 1000
dt = 0.5


# 1 . 1 ) Neuron / Network parameters (related to last slide of introduction)
n = 1000
pinh = 0.2
inh = (uniform(size=n)<pinh)
exc=logical_not(inh)
a = inh.choose(0.02,0.1)
b = 0.2 #RS , IB , FS: 0.2
c = -65 #RS , FS: −65 IB: −55
d = inh.choose(8,2)
tau_s = 10 # decay of synapses [ ms ]

#NEW recurrent parameters
width = pi / 4 # half−width of the orientation tuning
w = 0.0005 #average recurrent weight
pconn = 0.4 #recurrent connection prob
scaleEI = 2 # scaleI−>E
g_sc = 0.002 #scale of gamma
E = inh.choose(0,-85)
#NEW make weight matrix
W = zeros((n,n))
C = uniform(size=(n,n))
idx = nonzero(C<pconn) #sparse connectivity
W[idx] = gamma(w/g_sc,scale=g_sc,size=idx[0].size)
W[ix_(exc,inh)] *= scaleEI #submat indexing
theta = linspace(0,2*pi,n)
R = circulant(cos(theta))>cos(width)
W[:,exc] = where(R[:,exc],W[:,exc],0)
W = csr_matrix(W) #make row sparse

# 1 . 2 ) Input Parameters
tr= array([200,700] )/dt # stm time
rate_in = 2 # input rate
n_in = 100 # number of inputs
w_in = 0.07 # input weights
pconn_in = 0.1
C = uniform(size=(n,n_in))<pconn_in
W_in = C.choose(0,w_in)
W_in[int(n/2):,] = 0


# 2 ) reserve memory
T = int(ceil(tmax/dt))
v = zeros((T,n))
u = zeros((T,n))
v[0] = -70 # resting potential
u[0] = -14 # steady state
s_in = zeros( n_in ) # synaptic variable
E_in = zeros( n_in ) # rev potential
p_rate = dt * rate_in * 1e-3 # abbrev
s = zeros(n) #rec synapses


# 3 ) for−loop over time
for t in arange(T-1):
# 3 . 1 ) get input
    if t>tr[0] and t<tr[1] :
        # NEW: get input Poisson spikes
        p = uniform(size=n_in) < p_rate
    else :
        p = 0 # no input

    # NEW: calculate input current
    s_in = ( 1 - dt / tau_s ) * s_in + p
    I = dot( W_in , s_in  * E_in )
    I -= dot( W_in , s_in ) * v[t]

    fired = v[t]>=35

    # NEW: recurrent input
    s = (1 - dt / tau_s ) * s + fired
    Isyn = W.dot(s * E) - W.dot( s ) * v[t]
    I += Isyn # add to input vector

   
# 3 . 2 ) update ODE (ordinary differential equations)
    dv = (0.04 * v[t] + 5) * v[t] + 140 - u[t]
    v[t + 1] = v[t] + (dv + I) * dt
    du = a * (b * v [t] - u[t])
    u[t + 1] = u[t] + dt * du

# 3 . 3 ) spike!
    v[t][fired] = 35
    v[t + 1][fired] = c
    u[t + 1][fired] = u[t][fired] + d[fired]


# 4 ) plot voltage trace

tspk,nspk = nonzero(v==35)
idx_i = in1d(nspk,nonzero(inh)[0])
idx_e = logical_not(idx_i)

figure()
plot(tspk[idx_e]*dt,nspk[[idx_e]],'k.',label="Exc.",markersize=2)
plot(tspk[idx_i]*dt,nspk[[idx_i]],'r.',label="inh.",markersize=2)
xlabel("Time [ms]")
ylabel("Neuron Number [#]")
xlim(0,tmax)
title("An connected network of neurons")
legend(loc="upper right")
show()