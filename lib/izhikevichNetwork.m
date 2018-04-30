% function activity = izhikevichNetwork(T, n, distribution, connectivity) generates 
% spike trains from a network of Izhikevich's 2-D neuron model 
% 
% T: total simulation time (ms)
% n: number of neurons
% distribution: exponent on distribution of parameters c and d, which can bias the 
%   distribution toward adapting or bursting neurons (e.g. 2 biases toward adapting; 
%   0.1 biases toward bursting)
% connectivity: controls strength of excitatory connections (e.g. 0.5 is used by
%   Izhikevich, 0.8 produces more synchrony)
% 
% activity: a structure including: firings (an nX2 matrix with firing times in the 
%   first column and neuron #s in the second; spikes (spike times (s) with one
%   neuron per row); COV (coefficients of variation of each neuron's ISI)
% 
% function [activity, voltage] = izhikevichNetwork(..., vRange)
% 
% vRange: spread of random portion of initial voltage (default 0 causes
%   synchronous burst at startup). Voltages have range -65 to -65+vRange,
%   in mV
% 
% voltage: voltage history (asking for this slows things down considerably)
%
% Created by Eugene M. Izhikevich, February 25, 2003. 
%
% THIS CODE WAS ORIGINALLY TRANSCRIBED FROM IZHIKEVICH, 2003, IEEE TRANS
% NEURAL NETWORKS, WITH MINOR CHANGES BY BRYAN TRIPP, 2005 

function [activity, varargout] = izhikevichNetwork(T, n, distribution, connectivity, varargin)

    vrange = 0;
    if nargin > 4 % randomly initialize voltage 
        vrange = varargin{1};
    end
    
	% Excitatory neurons    Inhibitory neurons
	Ne = round(n*0.8);      Ni = round(n*0.2);
	re=rand(Ne,1);          ri=rand(Ni,1);
	a=[0.02*ones(Ne,1);     0.02+0.08*ri];
	b=[0.2*ones(Ne,1);      0.25-0.05*ri];
	c=[-65+15*re.^distribution;      -65*ones(Ni,1)];
	d=[8-6*re.^distribution;         2*ones(Ni,1)];
	S=[connectivity*rand(Ne+Ni,Ne),  -1*rand(Ne+Ni,Ni)];
	v=-65*ones(Ne+Ni,1) + vrange*rand(Ne+Ni,1); % Initial values of v
	u=b.*v; % Initial values of u
	firings=[]; % spike timings
    voltage = []; %voltage history
	for t=1:T 
        I=[5*randn(Ne,1);2*randn(Ni,1)]; % thalamic input
		fired=find(v>=30); % indices of spikes
		firings=[firings; t+0*fired,fired];
		v(fired)=c(fired);
		u(fired)=u(fired)+d(fired);
		I=I+sum(S(:,fired),2);
		v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
		v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
		u=u+a.*(b.*v-u);                 % stability
        
        if nargout > 1 %save voltage
            voltage = [voltage v];
        end
	end;
    voltage(find(voltage>30)) = 30;

    % put firing history in the form used elsewhere in this code
    spikes = zeros(Ne+Ni,floor(T/10));
    cov = zeros(1,Ne+Ni);
    for i = 1:Ne+Ni
        indices = find(firings(:,2) == i);
        ispikes = firings(indices,1)';
        spikes(i,1:length(indices)) = ispikes;
        isi = diff(ispikes);
        if length(isi) > 0
            cov(i) = std(isi) / mean(isi);
        end
    end
    
    activity = struct('firings', firings, 'spikes', spikes/1000, 'COV', cov);
    
    if nargout > 1
        varargout{1} = voltage;
    end
