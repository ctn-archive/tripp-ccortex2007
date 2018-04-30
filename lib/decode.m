% function [weights, err, t] = decode(signal, dt, spikes, randSD, bandSD, band, noiseRate, nt, ne, interval)
% 
% signal: target signal for ensemble to approximate
% dt: time step (s)
% spikes: a matrix of spike times with spike times for each neuron in a row
% randSD: see jitter.m
% bandSD: see jitter.m
% band: see jitter.m
% noiseRate: rate of extra Poisson-refractory spikes introduced as noise
% nt: number of "training" trials (repeated presentations used to find optimal weights)
% ne: number of evaluation trials (with novel noise, used to evaluate error)
% interval: sampling interval of signal and current used for finding weights (e.g. with dt=.0002 
%       and interval=5, weights would be found using 1000 Hz sampling rate)
% 
% weights: optimal weights for approximating given signal
% err: mean-squared error of approximation
% time: elapsed time to run this function  
% 
% function [..., corruptError, examples, estimates] = decode(..., plotEstimates)
% 
% plotEstimates: plots target signal and each approximation
% corruptError: error recalculated with slightly perturbed weights (to test
%    effect of numerical precision limitations)
% examples: firing patterns of a single neuron over different training and
%    evaluation trials
% estimates: optimal approximations of target for each evaluation trial 

% IMPLEMENTATION NOTE: moderate performance gains can be had by using raw 
% spike signals instead of PSCs, to estimate deconvolution of given signal
% instead of signal itself. However, too much fussing is then needed with 
% large signals / numbers of neurons, to avoid running out of memory. The 
% method below (convolving each spike train with PSC kernel individually) 
% turns out to be not much slower, and is more straightforward. 

function [weights, err, t, varargout] = decode(signal, dt, spikes, randSD, bandSD, band, noiseRate, nt, ne, interval, varargin)  
    tic
    compositeSignal = [];
    compositeCurrent = [];   
    psc = PSC(dt);
    
    T = dt*length(signal);

    examples = [];
    for i = 1:nt
        jittered = jitter(spikes, randSD, bandSD, band); 
        noised = addNoiseSpikes(jittered, T, dt, noiseRate);
        
        current = getCurrent(noised, dt, length(signal), psc); 
        compositeSignal = [compositeSignal, signal(1:interval:end)];
        compositeCurrent = [compositeCurrent, current(:,1:interval:end)];
        
        example = noised(1,:);
        examples(i,1:length(example)) = example;
    end
    
    pack

    weights = optimalDecoders(compositeSignal, compositeCurrent, 0);
    
    corruption = (-1 + 2*(rand(size(weights))>.5)) * eps;
    err = zeros(ne,1);
    
    if (nargin > 10 & varargin{1})
        time = dt*(1:length(signal));    
        figure, hold on, plot(time, signal, 'k--') 
    end
    estimates = [];
    for i = 1:ne
        jittered = jitter(spikes, randSD, bandSD, band);
        noised = addNoiseSpikes(jittered, T, dt, noiseRate);
        
        current = getCurrent(noised, dt, length(signal), psc);
        weighted = (weights * ones(size(signal))) .* current; 
        estimate = sum(weighted);
        
        if (nargin > 10 & varargin{1}) 
            plot(time, sum(weighted), 'k');
        end
        
        err(i) = mean( (sum(weighted) - signal).^2 );
        
        if (nargout > 3)
            corruptEstimate = ((weights+corruption) * ones(size(signal))) .* current; 
            corruptError(i) = mean( (sum(corruptEstimate) - signal).^2 );
        end

        example = noised(1,:);
        examples(nt+i,1:length(example)) = example;
        
        estimates = [estimates; estimate];
    end
    
    if (nargout > 3)
        varargout{1} = corruptError;
    end
    
    if (nargout > 4)
        varargout{2} = examples;
    end
    
    if (nargout > 5) 
        varargout{3} = estimates;
    end
    
    t = toc;
    
function noised = addNoiseSpikes(spikes, T, dt, noiseRate) 
    n = size(spikes,1);
    if (noiseRate == 0) 
        noiseSpikes = [];
    else 
        [noiseSpikes cov] = genUncorrelated(n, T, dt, noiseRate, [1 0 0]);
    end
    allSpikes = sort([spikes noiseSpikes], 2);
    
    % Now we may have ISIs smaller than the refractory time.  
    absRT = .001;
    for i = 1:n
        noised(i,:) = space(allSpikes(i,:), absRT);
    end
    
    