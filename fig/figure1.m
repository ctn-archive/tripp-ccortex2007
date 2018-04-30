% This scipt produces figure 1 of the article. 

rand('state',0)

dt = .0002;
signals = load('signals_figure1.mat');

tic
[activity, voltage] = izhikevichNetwork(1000, 1000, 2, 0.5, 10);
toc

%a target current based on smoothed, scaled mean firing rate
sumSpikes = sum(voltage > 0, 1);
sumSpikes = sumSpikes / mean(sumSpikes);
SD = 10; %smooth with 10ms SD Gaussian
time = (-5*SD):1:(5*SD);
gaussian = (1 / (SD*(2*pi)^.5)) * exp(-time.^2/(2*SD^2));
signalD = conv(sumSpikes, gaussian);
signalD = signalD((5*SD):end-(5*SD)-1);
signalD = interp(signalD, 5);

%plot raster 
figure(1), hold on, plot(activity.firings(:,1)/1000,activity.firings(:,2),'k.')
set(gca, 'Visible', 'off')
ylabel('Neuron #')
set(gca, 'XTick', [])

%plot example membrane potential
figure(2), hold on, plot(dt*[1:1000], voltage(101,:), 'k')
set(gca, 'Visible', 'off')
plot([.14 .14], [0 20], 'k')
text([.135], [25], '20 mV')
ylabel('Mem. potential (mV)')

%find optimal weights for each signal
%TODO: switch to fixed sign
%TODO: reinstate all training trials even without jitter (seems to effect
%results for some reason?)
% iterations = 1000;
% [weightsA, errorHistoryA] = learnedDecoders(signals.signalA, activity.spikes(1:500,:), zeros(500,1), iterations, 0, 0, activity.spikes(501:1000,:));
% [weightsC, errorHistoryC] = learnedDecoders(signals.signalC, activity.spikes(1:500,:), zeros(500,1), iterations, 0, 0, activity.spikes(501:1000,:));
% [weightsD, errorHistoryD] = learnedDecoders(signalD, activity.spikes(1:800,:), zeros(800,1), iterations, 0, 0, activity.spikes(801:1000,:));
% [weightsA, errA, t] = decode(signals.signalA, dt, activity.spikes, [0 0 0], [0 0 0], [0 1], 0, 1, 1, 1); t
% [weightsC, errC, t] = decode(signals.signalC, dt, activity.spikes, [0 0 0], [0 0 0], [0 1], 0, 1, 1, 1); t
% [weightsD, errD, t] = decode(signalD, dt, activity.spikes, [0 0 0], [0 0 0], [0 1], 0, 1, 1, 1); t
[weightsA, errA, t] = decode(signals.signalA, dt, activity.spikes(1:800,:), [0 0 0], [0 0 0], [0 1], 0, 1, 1, 1); t
[weightsC, errC, t] = decode(signals.signalC, dt, activity.spikes(1:800,:), [0 0 0], [0 0 0], [0 1], 0, 1, 1, 1); t
[weightsD, errD, t] = decode(signalD, dt, activity.spikes(1:800,:), [0 0 0], [0 0 0], [0 1], 0, 1, 1, 1); t

%get current at each synapse
psc = PSC(dt);
% TODO: switch
% current = getCurrent(activity.spikes, dt, length(signals.signalA), psc);
current = getCurrent(activity.spikes(1:800,:), dt, length(signals.signalA), psc);

%plot weighted currents as estimates of target patterns
figure(3), hold on
time = dt*(1:length(signals.signalA));
estimateA = (weightsA * ones(size(time))) .* current;
estimateC = (weightsC * ones(size(time))) .* current;
estimateD = (weightsD * ones(size(time))) .* current;
h = line(time, sum(estimateA)); set(h, 'LineWidth', [2]), set(h, 'Color', [.75 .75 .75]);
h = line(time, signals.signalA); set(h, 'LineWidth', [2]), set(h, 'Color', [0 0 0]);
h = line(time, 2.7+sum(estimateC)); set(h, 'LineWidth', [2]), set(h, 'Color', [.75 .75 .75]);
h = line(time, 2.7+signals.signalC); set(h, 'LineWidth', [2]), set(h, 'Color', [0 0 0]);
h = line(time, 4+sum(estimateD)); set(h, 'LineWidth', [2]), set(h, 'Color', [.75 .75 .75]);
h = line(time, 4+signalD); set(h, 'LineWidth', [2]), set(h, 'Color', [0 0 0]);
plot([.85 .95], [-2.4 -2.4], 'k')
plot([.85 .85], [-2.4 -1.4], 'k')
text(.78, -1.4, '1 nA')
text(.96, -2.4, '100 ms')
set(gca, 'Visible', 'off')
ylabel('Current (nA)')
xlabel('Time (s)')

% uncomment code below to find error with spike jitter
% clear current;
% clear estimateA;
% clear estimateC;
% clear estimateD;
% [weightsJA, errJA, t, ce, ex] = decode(signals.signalA, dt, activity.spikes, [.002 0 0], [0 0 0], [0 1], 0, 32, 5, 5); t

