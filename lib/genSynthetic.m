% function [spikes, cov] = genSynthetic(n, T, dt, correlatedRate, epochRate, 
% epochWidth, coarseCorrelatedRate, template, threshold, uncorrelatedRate, mix) 
% produces n synthetic spike trains.  
% 
% n: number of spike trains
% T: length of spike trains (seconds) 
% dt: temporal resolution (seconds)
% correlatedRate: the mean rate of correlated firing (spikes per second)
% epochRate: mean rate of correlation epochs (as arg c in genCorrelated.m)
% epochWidth: SD of spike distribution around epoch centre (as arg width in
%    genCorrelated.m)
% coarseCorrelatedRate: mean rate of coarsely correlated firing (see
%   genCoarseCorrelated.m)
% template: the signal arg of genCoarseCorrelated
% threshold: the threshold arg of genCoarseCorrelated
% uncorrelatedRate: the mean rate of uncorrelated firing (spikes per second)
% mix: uncorrelated interspike intervals are drawn from a weighted average of exponential, 
%   gaussian, and bimodel distributions.  Mix is a 3-vector that specifies these
%   weights.  
% 
% spikes: a matrix of n rows, each row giving the spike times of a separate 
%    spike train at least T seconds long
% cov: coefficients of variation of each spike train
% 
% [spikes, cov] = genSynthetic(..., options)
% 
% options: optional parameters for uncorrelated part of firing patterns
%    (see uncorrelated.m)

function [spikes, cov] = genSynthetic(n, T, dt, correlatedRate, epochRate, epochWidth, coarseCorrelatedRate, template, threshold, uncorrelatedRate, mix, varargin) 
    [uncorrSpikes cov] = genUncorrelated(n, T, dt, uncorrelatedRate, mix, varargin);
    corrSpikes = genCorrelated(n, T, dt, correlatedRate, epochRate, epochWidth);    
    [coarseSpikes, firingRate] = genCoarseCorrelated(n, dt, template, threshold, coarseCorrelatedRate);
    allSpikes = sort([uncorrSpikes corrSpikes coarseSpikes], 2);
    
    % now we may have ISIs smaller than the refractory time 
    absRT = .001; % 'absolute' refractory time
    for i = 1:n
        spikes(i,:) = space(allSpikes(i,:), absRT);
    end

    % find coefficients of variation
    for i = 1:n
        indices = find(spikes(i,:) > 0);
        isi = diff(spikes(i,indices));
        if length(isi) > 0
            cov(i) = std(isi) ./ mean(isi);
        else 
            cov(i) = 0;
        end
    end
    
    