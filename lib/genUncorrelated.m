% function [spikes, cov] = genUncorrelated(n, T, dt, rate, mix) produces n uncorrelated spike
% trains.  
% 
% n: number of spike trains
% T: length of spike trains (seconds) 
% dt: temporal resolution (seconds)
% rate: is the mean firing rate (spikes per second)
% mix: Interspike intervals are drawn from a weighted average of exponential, 
%   gaussian, and bimodel distributions.  Mix is a 3-vector that specifies these
%   weights.  
% 
% spikes: a matrix of n rows, each row giving the spike times of a separate 
%    spike train at least T seconds long
% cov: coefficients of variation of spike trains
% 
% [spikes, cov] = genUncorrelated(..., options)
% 
% options: a struct optionally containing any of the following fields: 
%    RT: refractory time for near-Poisson process (s)
%    SD: standard deviation of Gaussian ISI distribution for regularly firing neurons (s)
%    meanSD: SD of mean ISIs across population of regularly firing neurons (s)
%    SPB: mean spikes per burst in bursting distribution
%    COV: bursting coefficient of variation 
%    isSD: standard deviation of ISI within bursts (s)
%    bwCOV: coefficient of variation of inter-burst interval
%    periods: explicit list of mean ISIs per neuron

function [spikes, cov] = genUncorrelated(n, T, dt, rate, mix, varargin)  
    %defaults (see doUncorrelated for definitions) ... 
    RT = .002; 
    SD = .002;
    meanSD = .005; 	
    SPB = 5;
    COV = 2;
    inSD = .001;
    bwCOV = 1/3;
    
    if (nargin > 5)
        options = varargin{1};
        if isfield(options, 'RT'); RT = options.RT, end
        if isfield(options, 'SD'); SD = options.SD, end
        if isfield(options, 'meanSD'); meanSD = options.meanSD, end        
        if isfield(options, 'SPB'); SPB = options.SPB, end
        if isfield(options, 'COV'); COV = options.COV, end
        if isfield(options, 'inSD'); inSD = options.inSD, end
        if isfield(options, 'bwCOV'); bwCOV = options.bwCOV, end        
    end
    
    if rate > 0
        m = 1/rate + meanSD*randn(n,1);    
    end
    if (nargin > 5)
        if isfield(options, 'periods')
            m = options.periods;  
            if length(m) ~= n
                error(sprintf('List of mean firing periods is not the right size: %i', size(m)))
            end
        end
    end
    
    if (rate > 0)
        ISI = doUncorrelated(n, T, dt, rate, mix, RT, SD, m, SPB, COV, inSD, bwCOV);
        %make sure all trains cover time 0->T (since ISIs are random, no
        %  given number of ISIs guarantees this)
        while (min(sum(ISI,2)) < T) 
            ISI = [ISI, doUncorrelated(n, T/4, dt, rate, mix, RT, SD, m, SPB, COV, inSD, bwCOV)];
        end
	
        spikes = cumsum(ISI,2);
        cov = std(ISI, 0, 2) ./ mean(ISI, 2);
    else 
        spikes = [];
        cov = [];
    end

% args as above except ... 
% m: a list of n mean firing periods (mean ISIs) per neuron
function ISI = doUncorrelated(n, T, dt, rate, mix, RT, SD, m, SPB, COV, inSD, bwCOV)
    time = 0:dt:.3; % interspike interval time basis 
    if rate < 30 & rate > 0
        time = 0:dt:(10/rate);
    end
    
    %exponential distribution with refractory period for Poisson process 
    refractoryPDF = zeros(1,RT/dt);
    exponentialPDF = rate * exp(-rate*time);
    exponentialPDF = [refractoryPDF exponentialPDF(1:end-length(refractoryPDF))];    
    
    %bimodal ditribution for busting 
    s = COV / rate; %standard deviation of ISI derived from desired coefficient of variation 
    inMean = (1/rate) - (s^2/(SPB-0))^.5;
    bwMean = SPB/rate - (SPB-1)*inMean; %mean inter-burst interval such that overall firing rate is as specified above
    bwSD = bwCOV * bwMean;
    burstPDF = ((SPB-1)/SPB) * (1 / (inSD * (2*pi)^.5)) * exp(-(time-inMean).^2/(2*inSD^2)) + (1/SPB) * (1 / (bwSD * (2*pi)^.5)) * exp(-(time-bwMean).^2/(2*bwSD^2));
    
    nSpikes = ceil(T*rate); %number of ISIs to generate

    ISI = [];
    for i = 1:n
        %gaussian distribution for repetitive spiking (different for each train)
        if (SD == 0)
            repetitivePDF = zeros(size(time));
            index = round(m(i)/dt);
            if index <= length(time); repetitivePDF(index) = 1/dt; end
        else 
            repetitivePDF = (1 / (SD * (2*pi)^.5)) * exp(-(time-m(i)).^2/(2*SD^2));
        end

        ISI = [ISI; randPDF(1, nSpikes, mix(1)*exponentialPDF + mix(2)*repetitivePDF + mix(3)*burstPDF, dt)];
    end
    
