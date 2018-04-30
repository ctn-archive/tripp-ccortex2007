% This script is used to determine appropriate numbers of training trials
% for decode.m.  Error in generalization trials generally is reduced with
% greater numbers of training trials.  We want to use enough training
% trials so that the error approaches what it would be with a large number 
% of examples (as we assume are available in life), but we don't want to
% use too many, because the runtime increases with the number of trials, 
% and eventually we run into memory problems.  

clear all

x = load('signals_trainingTrials.mat');
signal = x.signal;

dt = .0002;
T = dt * length(signal);
spikeRate = 30;
[spikes, cov] = genUncorrelated(500, T, dt, spikeRate, [1 0 0]);

trialsCases = [1 2 4 8 16 32 45 64 128];
jitterCases = [0 .001 .01 0 0 .01];
noiseRateCases = spikeRate * [0 0 0 .5 2 2];

for i = 1:length(jitterCases)
    for j = 1:length(trialsCases)
        [weights, err, t(i,j)] = decode(signal, dt, spikes, [jitterCases(i) 0 0], [0 0 0], [0 1], noiseRateCases(i), trialsCases(j), 5, 5);  
        e(i,j) = mean(err);
        t(i,j)
    end
end

save 'data_trainingTrials.mat'
