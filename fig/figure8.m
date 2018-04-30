% This scipt produces figure 8 of the article from saved data.

load 'data_learning.mat'

figure
itNum = 0:iterations;
loglog(itNum, error, 'b')
set(gca, 'NextPlot', 'add')
plot(itNum, errorJitter, 'c')
plot(itNum, errorFilter50, 'r')
plot(itNum, errorFilter500, 'k')
plot(itNum, errorFilterJitter, 'k')
set(gcf, 'Position', [360 669 308 265])

psc = PSC(dt);
components = getCurrent(spikes, dt, 1500, psc);
time = dt:dt:.3;

weighted = weights * ones(size(signal)) .* components;
figure, hold on, plot(time, sum(weighted), 'k'), plot(time, signal, 'k'), set(gcf, 'Position', [586 311 206 249])
set(gca, 'XLim', [0 .3]), set(gca, 'YLim', [0 1.5])

weighted = weightsFilter50 * ones(size(signal)) .* components;
figure, hold on, plot(time, sum(weighted), 'k'), plot(time, signal, 'k'), set(gcf, 'Position', [586 311 206 249])
set(gca, 'XLim', [0 .3]), set(gca, 'YLim', [0 1.5])

weighted = weightsJitter * ones(size(signal)) .* components;
figure, hold on, plot(time, sum(weighted), 'k'), plot(time, signal, 'k'), set(gcf, 'Position', [586 311 206 249])
set(gca, 'XLim', [0 .3]), set(gca, 'YLim', [0 1.5])
plot([.1 .2], [.2 .2], 'k')
plot([.1 .1], [.2 .7], 'k')

