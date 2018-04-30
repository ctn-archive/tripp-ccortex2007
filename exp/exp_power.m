% Finds the number of trials needed to find the minimal spike rate increases
% that would be needed to support reliable spiking. Assumptions as in article text.  

clear all

lambda = 0.1; %baseline spike rate per bin
bins = 200;
N = 10000; %total number of converging neurons 
n = [500 1000]; %number of converging neurons that contribute to each spike 
                % ... don't want < 500 because then we begin to have
                % neurons that don't contribute to any output spike                

ne = lambda * bins * n / N;

mistimeRate = [.001:.005:.09]; %error rate of post-synaptic neuron (per bin)
probExtraSpike = mistimeRate;
probMissingSpike = mistimeRate / lambda; %this is the probability of any intended spike getting dropped
z_baseline = norminv(1-probExtraSpike, 0, 1)
z_elevated = norminv(1-probMissingSpike, 0, 1)

presynapticRate = bins * lambda;

minimalElevatedRate = [];
for i = 1:length(n) 
    for j = 1:length(mistimeRate)
        threshold = lambda + z_baseline(j) * (lambda/N)^.5; % note the N (all neurons may combine to spike inappropriately)
        elevated = elevatedRate(lambda, n(i), N, threshold, probMissingSpike(j));
        minimalElevatedRate(i,j) = elevated;
        thresholds(i,j) = threshold;
        
        % find presynaptic error rate (to allow comparison with post-synaptic error rate) 
        presynapticExtraSpikeRate(i,j) = (bins-presynapticRate) * (1-pdf('poiss', 0, lambda)); %i.e. #supposedly non-spiking bins X prob >0 spikes/bin
        presynapticMissingSpikeRate(i,j) = ne(i) * pdf('poiss', 0, elevated) + (presynapticRate - ne(i)) * pdf('poiss', 0, lambda);
        presynapticErrorRate(i,j) = (presynapticExtraSpikeRate(i,j) + presynapticMissingSpikeRate(i,j)) / 2;
    end
end

% display some intermediate results
thresholds
minimalElevatedRate
presynapticExtraSpikeRate
presynapticMissingSpikeRate
presynapticErrorRate

rateDifference = minimalElevatedRate - lambda;

sigma = lambda ^ .5;

% Differences would be best detected with ANOVA, so we find # trials needed for sufficient ANOVA power, with effect sizes as calculated above 
% need #elevated trials, lambda, lamba_elevated, bins
power = .8; %desired power 
for i = 1:length(n)
    for j = 1:length(mistimeRate)
        tic 
        trials(i,j) = anovaTrials(bins, lambda, minimalElevatedRate(i,j), ne(i), power)
        toc
        save 'data_power.mat'
    end    
end

