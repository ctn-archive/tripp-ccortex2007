% This script experiments with sensitivity to ISI jitter with varying 
% population sizes and firing rates. 

clear all

jitter = [0 .001 .002 .004 .008 .016];

% population size, mean firing rate, Poisson mix coefficient, regular coefficient, burst coefficient, Poisson refractory time
cases = [ ... 
        500 2 1 0 0 .002; ...
        500 5 1 0 0 .002; ...
        500 10 1 0 0 .002; ...
        500 20 1 0 0 .002; ... 
        500 40 1 0 0 .002; ...
        500 80 1 0 0 .002; ...
        500 160 1 0 0 .002; ...
        500 320 1 0 0 .002; ...        
        500 2 1 0 0 0; ...
        500 5 1 0 0 0; ...
        500 10 1 0 0 0; ...
        500 20 1 0 0 0; ... 
        500 40 1 0 0 0; ...
        500 80 1 0 0 0; ...
        500 160 1 0 0 0; ...
        500 320 1 0 0 0; ... 
        250 40 1 0 0 .002; ...
        750 40 1 0 0 .002; ...
        1000 40 1 0 0 .002; ...
        1250 40 1 0 0 .002; ...
        1500 40 1 0 0 .002; ...
        ]
    
dt = .0002;

T = .3; 
covT = 20; % longer simulations used for COV calculations 

x = load('signals_figure4.mat');
signals = x.signals;

meanCOV = [];
meanE = [];
sdE = [];

for i = 1:size(cases,1)
    params = struct('RT', cases(i,6));
    
    % generate long firing patterns to get COV estimate
    [spikes, cov] = genUncorrelated(20, covT, dt, cases(i,2), cases(i,3:5), params);
    meanCOV = [meanCOV; mean(cov)]

    % We have to use a 1ms intervals with large populations, to avoid running out of memory. However with 
    % 1ms intervals, firing rates < 10Hz, and no jitter, about one case in 10 has abnormally high error (about 
    % an order of magnitude greater than usual). This doesn't happen with .2ms intervals. We assume this is a 
    % numerical problem, and try to avoid this case. 
    if cases(i,1) > 600
        interval = 5;
    else 
        interval = 1;
    end
    
    for j = 1:length(jitter) 
        for k = 1:size(signals,1)            
            % generate test firing patterns 
            [spikes, cov] = genUncorrelated(cases(i,1), T, dt, cases(i,2), cases(i,3:5), params);
            
            % find weights to approximate target current signal
            [weights, err, t] = decode(signals(k,:), .0002, spikes, [jitter(j) 0 0], [0 0 0], [0 10], 0, 32, 5, interval); 
            t
            meanE(i,j,k) = mean(err)
            sdE(i,j,k) = std(err);  
        end
    end
    save 'data_populationSizeRate.mat'
end

