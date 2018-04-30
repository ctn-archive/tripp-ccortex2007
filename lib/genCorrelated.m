% function spikes = genCorrelated(n, T, dt, rate) produces a set of spike trains
% that are correlated around Poisson-distributed correlation times
% 
% n: number of spike trains
% T: length of spike trains (seconds)
% dt: time step (seconds)
% rate: mean spike rate (spikes per second; maximum 75
% 
% Returns a matrix of spike times where each row corresponds to a spike train
% 
% function spikes = genCorrelated(n, T, dt, rate, c, width, skewness)  
% 
% c: number of points in time at which trains are correlated (default 75) 
% width: SD of distribution of spikes around these points
% skewness: of distributions of spikes around correlation points (1Xn array)
% 
% Note: this method is adapted from Benucci, Verschure, & Konig, 2004 

function spikes = genCorrelated(n, T, dt, rate, varargin)
    epochRate = 75;
    sd = .002;
    skewness = rand(1,n) * 6 - 3;
    if (nargin > 4); epochRate = varargin{1}; end
    if (nargin > 5); sd = varargin{2}; end
    if (nargin > 6); skewness = varargin{3}; end       
    
    if (rate > epochRate) 
        error('Firing rate can not be higher than epoch rate')
    end
    
    interEpochIntervalPDF = epochRate * exp(-epochRate*[0:dt:15/epochRate]);
    epochIntervals = randPDF(1, 2*T*epochRate, interEpochIntervalPDF, dt);
    epochTimes = cumsum(epochIntervals);
    epochTimes = epochTimes(find(epochTimes <= T));
    
    probPerEpoch = rate / epochRate;
    spikes = zeros(n,length(epochTimes));
    for i = 1:n
        ithNeuronSpikes = [];
        basis = [(-sd*4):dt:(sd*4)];
        normalPDF = (1 / (sd * (2*pi)^.5)) * exp(-(basis).^2/(2*sd^2));
        skew = (1/2) * (1 + erf(skewness(i)*basis/(sd*2^.5)));
        skewedPDF = 2*normalPDF.*skew;
        shifts = randPDF(1,length(epochTimes),skewedPDF,dt) + (basis(1) - dt); %latter term centres on zero
        for e = 1:length(epochTimes)
            if (probPerEpoch >= rand)
                ithNeuronSpikes = [ithNeuronSpikes, epochTimes(e) + shifts(e)];
            end
        end
        spikes(i,1:length(ithNeuronSpikes)) = ithNeuronSpikes;
    end
