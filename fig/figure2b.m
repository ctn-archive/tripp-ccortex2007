% This script produces figure 2 panel E of the article 

frequencies = [.125 .25 .5 1 2 4 8 16 32 64 128 256 512];

mixes = [0 1 0; ...
       1 15 0; ...
       1 6 0; ...
       1 3 0; ...
       1 1 0; ...
       2 1 0; ...
       1 0 0; ...
       1 0 .1; ...
       1 0 .25; ...
       1 0 .6; ...
       1 0 1.5; ...
       1 0 5];
   
dt = .0002;
nStats = 20; % # neurons used to estimate stats 
n = 500; % # neurons used for experiment
T = .5;
       
COVs = zeros(size(mixes,1), 1);
errors = zeros(size(mixes,1), length(frequencies));

for i = 1:size(mixes,1)
    mix = mixes(i,:);
    [spikes, COV] = genUncorrelated(nStats, 10, dt, 30, mix, struct('SD', .0025, 'meanSD', .0025)); 
    COVs(i) = mean(COV)

    [spikes, COV] = genUncorrelated(n, T+.1, dt, 30, mix, struct('SD', .0025, 'meanSD', .0025));
    spikes = spikes - .1; %to avoid dead space at beginning with periodic spiking
    
    phaseErrors = frequency(spikes, frequencies, 0);
    errors(i,:) = mean(phaseErrors, 2)'
end

figure, hold on
mesh(frequencies, COVs, errors);
set(gca, 'XScale', 'log')
set(gca, 'ZScale', 'log')

xlabel('Frequency (Hz)')
ylabel('COV')
zlabel('MSE')

set(gcf, 'Position', [360 669 308 265])
