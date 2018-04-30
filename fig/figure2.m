% This script produces figure 2 of the article 

clear all

rand('state',0)

signal = load('signals_figure2.mat')

dt = .0002;
nStats = 20; % # neurons used to estimate stats 
n = 500; % # neurons used for experiment
time = dt*[1:length(signal.signal)]; 
T = max(time); 

% generate long patterns to estimate coefficients of variations
[spikes, covA] = genUncorrelated(nStats, 10, dt, 30, [0 1 0], struct('SD', .0025, 'meanSD', .0025)); mean(covA)
[spikes, covB] = genUncorrelated(nStats, 10, dt, 30, [1 0 0]); mean(covB)
[spikes, covC] = genUncorrelated(nStats, 10, dt, 30, [1 0 5]); mean(covC)

% generate larger populations (same parameters) for experiment
[spikesA, cov] = genUncorrelated(n, 2*T, dt, 30, [0 1 0], struct('SD', .0025, 'meanSD', .0025));
[spikesB, cov] = genUncorrelated(n, 2*T, dt, 30, [1 0 0]);
[spikesC, cov] = genUncorrelated(n, 2*T, dt, 30, [1 0 5]);

% generate ISIs for population with flat firing rate distribution 
ratesD = 60 * (1:n) / n;
periodsD = 1 ./ ratesD;
periodsD = periodsD(randperm(n));
startsD = rand(size(periodsD)) .* periodsD;
for i = 1:n
    spikes = startsD(i);
    while spikes(end) < 2*T 
        spikes = [spikes spikes(end) + periodsD(i)];
    end
    spikesD(i,1:length(spikes)) = spikes;
end

% find optimal synaptic weights
[weightsA, errA, t] = decode(signal.signal, dt, spikesA, [0 0 0], [0 0 0], [0 1], 0, 1, 1, 1); t, errA
[weightsB, errB, t] = decode(signal.signal, dt, spikesB, [0 0 0], [0 0 0], [0 1], 0, 1, 1, 1); t, errB
[weightsC, errC, t] = decode(signal.signal, dt, spikesC, [0 0 0], [0 0 0], [0 1], 0, 1, 1, 1); t, errC
[weightsD, errD, t] = decode(signal.signal, dt, spikesD, [0 0 0], [0 0 0], [0 1], 0, 1, 1, 1); t, errD

psc = PSC(dt);

% get current at each synapse
currentA = getCurrent(spikesA, dt, 2*length(time), psc);
currentB = getCurrent(spikesB, dt, 2*length(time), psc);
currentC = getCurrent(spikesC, dt, 2*length(time), psc);
currentD = getCurrent(spikesD, dt, 2*length(time), psc);

% calculate estimate as weighted sum of post-synaptic currents 
estimateA = (weightsA * ones(size(time))) .* currentA(:,1:length(time));
estimateB = (weightsB * ones(size(time))) .* currentB(:,1:length(time));
estimateC = (weightsC * ones(size(time))) .* currentC(:,1:length(time));
estimateD = (weightsD * ones(size(time))) .* currentD(:,1:length(time));

% plot example rasters and signal estimates 
figure(1), hold on
h = line(time, sum(estimateA)); set(h, 'LineWidth', [2]), set(h, 'Color', [.75 .75 .75]);
h = line(time, signal.signal); set(h, 'LineWidth', [2]), set(h, 'Color', [0 0 0]);
for i = 1:20
    indices = find(spikesA(i,:) > 0 & spikesA(i,:) < T);
    plot(spikesA(i,indices), -2.5+(i/5)*ones(size(indices)), 'k.');
end
set(gcf, 'Position', [586 311 206 249])
set(gca, 'XLim', [0 .5])

figure(3), hold on
h = line(time, sum(estimateB)); set(h, 'LineWidth', [2]), set(h, 'Color', [.75 .75 .75]);
h = line(time, signal.signal); set(h, 'LineWidth', [2]), set(h, 'Color', [0 0 0]);
for i = 1:20
    indices = find(spikesB(i,:) > 0 & spikesB(i,:) < T);
    plot(spikesB(i,indices), -2.5+(i/5)*ones(size(indices)), 'k.');
end
set(gcf, 'Position', [586 311 206 249])
set(gca, 'YTick', [])
set(gca, 'XLim', [0 .5])

figure(4), hold on
h = line(time, sum(estimateC)); set(h, 'LineWidth', [2]), set(h, 'Color', [.75 .75 .75]);
h = line(time, signal.signal); set(h, 'LineWidth', [2]), set(h, 'Color', [0 0 0]);
for i = 1:20
    indices = find(spikesC(i,:) > 0 & spikesC(i,:) < T);
    plot(spikesC(i,indices), -2.5+(i/5)*ones(size(indices)), 'k.');
end
plot([.53 .63], [0 0], 'k') % scale bars
plot([.53 .53], [0 1], 'k')
set(gcf, 'Position', [840 471 266 268])
set(gca, 'YTick', [])
set(gca, 'XLim', [0 .65])

figure(5), hold on
h = line(time, sum(estimateD)); set(h, 'LineWidth', [2]), set(h, 'Color', [.75 .75 .75]);
h = line(time, signal.signal); set(h, 'LineWidth', [2]), set(h, 'Color', [0 0 0]);
for i = 1:20
    indices = find(spikesD(i,:) > 0 & spikesD(i,:) < T);
    plot(spikesD(i,indices), -2.5+(i/5)*ones(size(indices)), 'k.');
end
set(gcf, 'Position', [586 311 206 249])
set(gca, 'YTick', [])
set(gca, 'XLim', [0 .5])

% find frequency content of error for first example
error = sum(estimateA) - signal.signal;
f = fft(error);
p = f(1:1250);
ap = p .* conj(p);
figure(2), plot([2:2:2500], cumsum(ap))
max(cumsum(ap)) / 4

% plot frequency content of first 5 principle components 
npc = 5; 
lenfft = 5000; % not too long for PCA, 1Hz increments, assuming dt=.0002;

[PC, score, latent, tsquare] = princomp(currentA); %gives components in columns
PCS = fft(PC(1:lenfft,1:npc));
PCS = PCS .* conj(PCS);
figure(6), plot(0:100, PCS(1:101,:))
xlabel('Frequency (Hz)')
set(gcf, 'Position', [286 311 206 249])
set(gca, 'YTick', [])

[PC, score, latent, tsquare] = princomp(currentB); 
PCS = fft(PC(1:lenfft,1:npc));
PCS = PCS .* conj(PCS);
figure(7), plot(0:100, PCS(1:101,:))
xlabel('Frequency (Hz)')
set(gcf, 'Position', [286 311 206 249])
set(gca, 'YTick', [])

[PC, score, latent, tsquare] = princomp(currentC); 
PCS = fft(PC(1:lenfft,1:npc));
PCS = PCS .* conj(PCS);
figure(8), plot(0:100, PCS(1:101,:))
xlabel('Frequency (Hz)')
set(gcf, 'Position', [286 311 206 249])
set(gca, 'YTick', [])

[PC, score, latent, tsquare] = princomp(currentD); 
PCS = fft(PC(1:lenfft,1:npc));
PCS = PCS .* conj(PCS);
figure(9), plot(0:100, PCS(1:101,:))
xlabel('Frequency (Hz)')
set(gcf, 'Position', [286 311 206 249])
set(gca, 'YTick', [])

% % plot error as a function of frequency for cases A-D
% frequencies = [.125 .25 .5 1 2 4 8 16 32 64 128 256 512];
% freqErrorA = frequency(spikesA, frequencies, 0);
% freqErrorB = frequency(spikesB, frequencies, 0);
% freqErrorC = frequency(spikesC, frequencies, 0);
% freqErrorD = frequency(spikesD, frequencies, 0);
% figure(10)
% loglog(frequencies, mean(freqErrorA, 2), 'r');
% set(gca, 'NextPlot', 'add');
% loglog(frequencies, mean(freqErrorB, 2), 'm');
% loglog(frequencies, mean(freqErrorC, 2), 'b');
% loglog(frequencies, mean(freqErrorD, 2), 'g');
% for i = 1:length(frequencies)
%     h = line([frequencies(i) frequencies(i)], [mean(freqErrorA(i,:))-std(freqErrorA(i,:)) mean(freqErrorA(i,:))+std(freqErrorA(i,:))]); 
%     set(h, 'Color', 'r'), set(h, 'LineWidth', 2)
%     h = line([frequencies(i) frequencies(i)], [mean(freqErrorB(i,:))-std(freqErrorB(i,:)) mean(freqErrorB(i,:))+std(freqErrorB(i,:))]);
%     set(h, 'Color', 'm'), set(h, 'LineWidth', 2)
%     h = line([frequencies(i) frequencies(i)], [mean(freqErrorC(i,:))-std(freqErrorC(i,:)) mean(freqErrorC(i,:))+std(freqErrorC(i,:))]);
%     set(h, 'Color', 'b'), set(h, 'LineWidth', 2)
%     h = line([frequencies(i) frequencies(i)], [mean(freqErrorD(i,:))-std(freqErrorD(i,:)) mean(freqErrorD(i,:))+std(freqErrorD(i,:))]);
%     set(h, 'Color', 'g'), set(h, 'LineWidth', 2)
% end
% get(h)
% xlabel('Frequency (Hz)')
% ylabel('MSE')
% set(gcf, 'Position', [360 669 308 265])
% set(gca, 'XLim', [.05 1000])
% set(gca, 'YLim', [0.0005 1])


% Note: the remaining panels of figure 2 can be produced with figure2b.m 
% and see figure2c.m 
