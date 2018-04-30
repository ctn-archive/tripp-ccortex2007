% this script produces figure 6 of the article

% plot error vs correlation summary 
jitterIndex = 4;
figure(1), plotCorr('data_correlation.mat', jitterIndex)
set(gcf, 'Position', [360 669 308 265])

% plot examples ... 
clear all

rand('state',0)

signal = load('signals_figure2.mat')

dt = .0002;
nStats = 20;
n = 500;
time = dt*[1:length(signal.signal)];
T = max(time);
corrT = 20;

% a few long spike trains for calculating statistics ... 
%poisson correlated 
[spikes cov] = genSynthetic(nStats, corrT, dt, 40, 100, .003, 0, sin([dt:dt:corrT]*70), 0, 0, [1 0 0]);
xc = ensembleCrossCorr(spikes, .001, corrT);
sprintf('Correlation: %2.2f +/- %2.2f   COV: %2.2f +/- %2.2f', mean(xc(:,1)), std(xc(:,1)), mean(cov), std(cov))

%alpha correlated 
[spikes cov] = genSynthetic(nStats, corrT, dt, 0, 75, .003, 40, sin([dt:dt:corrT]*2*pi*10), -1, 0, [1 0 0]);
xc = ensembleCrossCorr(spikes, .001, corrT);
sprintf('Correlation: %2.2f +/- %2.2f   COV: %2.2f +/- %2.2f', mean(xc(:,1)), std(xc(:,1)), mean(cov), std(cov))

%gamma correlated
[spikes cov] = genSynthetic(nStats, corrT, dt, 0, 75, .003, 40, sin([dt:dt:corrT]*2*pi*55), -1, 0, [1 0 0]);
xc = ensembleCrossCorr(spikes, .001, corrT);
sprintf('Correlation: %2.2f +/- %2.2f   COV: %2.2f +/- %2.2f', mean(xc(:,1)), std(xc(:,1)), mean(cov), std(cov))

%many shorter spike trains (same parameters) for simulations ... 
[spikesA, cov] = genSynthetic(n, T, dt, 40, 100, .003, 0, sin([dt:dt:T]*70), 0, 0, [1 0 0]);
[spikesB, cov] = genSynthetic(n, T, dt, 0, 75, .003, 40, sin([dt:dt:T]*2*pi*10), -1, 0, [1 0 0]);
[spikesC, cov] = genSynthetic(n, T, dt, 0, 75, .003, 40, sin([dt:dt:T]*2*pi*55), -1, 0, [1 0 0]);

%find optimal synaptic weights 
[weightsA, errA, t] = decode(signal.signal, dt, spikesA, [0 0 0], [0 0 0], [0 1], 0, 1, 1, 1); t, errA
[weightsB, errB, t] = decode(signal.signal, dt, spikesB, [0 0 0], [0 0 0], [0 1], 0, 1, 1, 1); t, errB
[weightsC, errC, t] = decode(signal.signal, dt, spikesC, [0 0 0], [0 0 0], [0 1], 0, 1, 1, 1); t, errC

psc = PSC(dt);

%get current at each synapse
currentA = getCurrent(spikesA, dt, length(time), psc);
currentB = getCurrent(spikesB, dt, length(time), psc);
currentC = getCurrent(spikesC, dt, length(time), psc);

%target pattern approximations as weighted sums of synaptic currents 
estimateA = (weightsA * ones(size(time))) .* currentA;
estimateB = (weightsB * ones(size(time))) .* currentB;
estimateC = (weightsC * ones(size(time))) .* currentC;

%plot rasters and estimates
nRasters = 20;
figure(2)
h11 = axes, hold on, set(gca, 'Visible', 'off'), set(gca, 'XLim', [0 .5])
h = line(time, sum(estimateA)); set(h, 'LineWidth', [2]), set(h, 'Color', [.75 .75 .75]);
h = line(time, signal.signal); set(h, 'LineWidth', [2]), set(h, 'Color', [0 0 0]);
h12 = axes, hold on, set(gca, 'Visible', 'off'), set(gca, 'XLim', [0 .5])
for i = 1:nRasters
    indices = find(spikesA(i,:) > 0 & spikesA(i,:) < T);
    plot(spikesA(i,indices), i*ones(size(indices)), 'k.');
end
set(gcf, 'Position', [586 311 206 249])

figure(3)
h21 = axes, hold on, set(gca, 'Visible', 'off'), set(gca, 'XLim', [0 .5])
h = line(time, sum(estimateB)); set(h, 'LineWidth', [2]), set(h, 'Color', [.75 .75 .75]);
h = line(time, signal.signal); set(h, 'LineWidth', [2]), set(h, 'Color', [0 0 0]);
h22 = axes, hold on, set(gca, 'Visible', 'off'), set(gca, 'XLim', [0 .5])
for i = 1:nRasters
    indices = find(spikesB(i,:) > 0 & spikesB(i,:) < T);
    plot(spikesB(i,indices), i*ones(size(indices)), 'k.');
end
set(gcf, 'Position', [586 311 206 249])

figure(4)
h31 = axes, hold on, set(gca, 'Visible', 'off'), set(gca, 'XLim', [0 .65])
h = line(time, sum(estimateC)); set(h, 'LineWidth', [2]), set(h, 'Color', [.75 .75 .75]);
h = line(time, signal.signal); set(h, 'LineWidth', [2]), set(h, 'Color', [0 0 0]);
h32 = axes, hold on, set(gca, 'Visible', 'off'), set(gca, 'XLim', [0 .65])
for i = 1:nRasters
    indices = find(spikesC(i,:) > 0 & spikesC(i,:) < T);
    plot(spikesC(i,indices), i*ones(size(indices)), 'k.');
end
axes(h31)
plot([.53 .63], [0 0], 'k')
plot([.53 .53], [0 1], 'k')
text(.55, -1, '100 ms')
text(.6, 1, '1 nA')
set(gcf, 'Position', [581 489 265 249])
