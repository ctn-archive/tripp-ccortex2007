% This script experiments with learning

clear all

dt = .0002;
x = load('signals_figure4.mat');
signal = x.signals(2,:);

[spikes, cov] = genUncorrelated(500, .3, dt, 30, [1 0 0]);
[statspikes, cov] = genUncorrelated(20, 10, dt, 200, [0 1 0], struct('SD', .001, 'meanSD', .002));
mean(cov)
std(cov)

iterations = 1000;

[weights, error] = learnedDecoders(signal, spikes, zeros(500,1), iterations, 0, 0);
[weightsJitter, errorJitter] = learnedDecoders(signal, spikes, zeros(500,1), iterations, 0.004, 0);
[weightsFilter50, errorFilter50] = learnedDecoders(signal, spikes, zeros(500,1), iterations, 0, 0.05);
[weightsFilter500, errorFilter500] = learnedDecoders(signal, spikes, zeros(500,1), iterations, 0, 0.5);
 
[weightsFilterJitter, errorFilterJitter] = learnedDecoders(signal, spikes, zeros(500,1), iterations, 0.004, 0.05);
save 'data_learning.mat';
