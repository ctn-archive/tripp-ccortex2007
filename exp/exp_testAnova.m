% Script to test that type-I error is roughly maintained with Poisson 
% variables instead of normal variables (within ranges of interest)

clear all 

bins = 200;
lambda = .1;
ne = 0;
alpha = .05;

trials = [10 31 100 315 1000 3150]; %make sure it works over a wide range of #s of trials per experiment

for i = 1:length(trials)
    for j = 1:5
        power(i,j) = anovaPowerExperiment(bins, lambda, lambda, 0, trials(i), alpha)
    end
end

p = mean(power');
sdp = std(power');

figure
semilogx(trials, p, 'k.')
set(gcf, 'Position', [360 669 308 265])
set(gca, 'NextPlot', 'add')

for i = 1:length(trials)
    semilogx([trials(i) trials(i)], [p(i)-sdp(i) p(i)+sdp(i)], 'k')
end

set(gca, 'YLim', [0 .1])
set(gca, 'XLim', [min(trials)/2 max(trials)*2])
