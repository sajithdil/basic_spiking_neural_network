% spnet.m: Spiking network with axonal conduction delays and STDP
% Created by Eugene M.Izhikevich. February 3, 2004
M=100; % number of synapses per neuron
D=20; % maximal conduction delay
% excitatory neurons % inhibitory neurons % total number
Ne=800; Ni=200; N=Ne+Ni; % RS -> a =0.02/d=8, FS-> a = 0.1/d=2 , 
a=[0.02*ones(Ne,1); 0.1*ones(Ni,1)];% Excitatory Neurons = cortical pyramidal neurons RS firing / inhibitory cortical interneurons fast spiking
d=[ 8*ones(Ne,1); 2*ones(Ni,1)];
sm=10; % maximal synaptic strength . so that it would not grow beyond
post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]);% ceil fn rounds off values to the highest poitive integer/ post-synaptic targets
s=[6*ones(Ne,M);-5*ones(Ni,M)]; % initial synaptic weights
sd=zeros(N,M); % synaptic weight derivative initiation
for i=1:N
    if i<=Ne
        for j=1:D
            delays{i,j}=M/D*(j-1)+(1:M/D); % Neuron i is having a delay of j / excitatory neurons have 1 to D=20 milisecond delays
        end;
    else
        delays{i,1}=1:M;% 1 ms delay to all inhibitory neurons
    end;
    pre{i}=find(post==i&s>0); % amma gahai bat eken, find the indices macthing 1 to 1000. One at a time \  all excitatory neurons presynaptic to a given neuron
    aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);% auxilary table of indices needed to speed up STDP implementation
end;
STDP = zeros(N,1001+D);% 1000 by 1021 matrix of zeros
v = -65*ones(N,1); % initial values izhi neuron model
u = 0.2.*v; % initial values izhi neuron model
firings=[-D 0]; % spike timings
for sec=1:60*60*24 % simulation of 1 day
    for t=1:1000 % simulation of 1 sec (since all values are in milisecond scale)
        I=zeros(N,1); % Set of values to excite the input neurons
        I(ceil(N*rand))=20; % random thalamic input of 20mV (spikes one random neuron each iteration)
        fired = find(v>=30); % indices of fired neurons(Find the neurons that have fired)
        v(fired)=-65; % resetting the fired neuron | typical c = -65 according to izhi model
        u(fired)=u(fired)+d(fired); % u = memebrain recovery variable | resetting formula
        STDP(fired,t+D)=0.1;% each time a neuron fires, stdp is reset to 0.1
        for k=1:length(fired) % for the length of the indecies fired
            sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});%
        end;
        firings=[firings;t*ones(length(fired),1),fired];
        k=size(firings,1);
        while firings(k,1)>t-D
            del=delays{firings(k,2),t-firings(k,1)+1};
            ind = post(firings(k,2),del); 
            I(ind)=I(ind)+s(firings(k,2), del)';
            sd(firings(k,2),del)=sd(firings(k,2),del)-1.2*STDP(ind,t+D)';
            k=k-1;
        end;
        v=v+0.5*((0.04*v+5).*v+140-u+I); % for numerical
        v=v+0.5*((0.04*v+5).*v+140-u+I); % stability time
        u=u+a.*(0.2*v-u); % step is 0.5 ms
        STDP(:,t+D+1)=0.95*STDP(:,t+D); % tau = 20 ms
    end;
    plot(firings(:,1),firings(:,2),'.');
    axis([0 1000 0 N]); drawnow;
    STDP(:,1:D+1)=STDP(:,1001:1001+D);
    ind = find(firings(:,1) > 1001-D);
    firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
    s(1:Ne,:)=max(0,min(sm,0.01+s(1:Ne,:)+sd(1:Ne,:)));
    sd=0.9*sd;
end;