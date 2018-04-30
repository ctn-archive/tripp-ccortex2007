% This script produces figure 3 of the article

clear all

dt = .0002;
T = 0.5;
x = load('signals_figure3.mat');
signal = x.signal3;

rand('state',1)

% generate base pattern
[spikes, cov] = genUncorrelated(1500, T, dt, 20, [1 0 0]);

nt = 32;
ne = 2;

% find optimal weights under noise
tic
[weights, err, t, corrupt, examples, estimates] = decode(signal, dt, spikes, [.02 0 0], [0 0 0], [0 1], 20, nt, ne, 5, 1);
toc

% plot estimate
plot([.05 .05], [9 11], 'k')
xlabel('Time (s)')
ylabel('Current (nA)')

% plot rasters
figure, hold on
for i = 1:nt
    indices = find(examples(i,:) > 0 & examples(i,:) < T);    
    plot(examples(i,indices), i*ones(size(indices)), 'k.');
end
for i = (nt+1):(nt+ne)
    indices = find(examples(i,:) > 0 & examples(i,:) < T); 
    plot(examples(i,indices), i*ones(size(indices)), 'ko');    
end
xlabel('Time (s)')
ylabel('Trial #')

% plot histogram 
binSize = .01;
histogram = zeros(1,T/binSize);
for i = 1:size(examples,1)
    for j = 1:size(examples,2)
        bin = ceil(examples(i,j) / binSize);
        if (bin > 0 & bin < length(histogram))
            histogram(bin) = histogram(bin) + 1;
        end
    end
end
figure, bar((1/binSize) * histogram / size(examples,1))
ylabel('Firing rate (Hz)')

%feed into Hodgkin-Huxley model
global HH_T HH_I
HH_T = 0:.2:500;
HH_I = estimates(1,:);
[T1,r1] = ode45('hh', [0 499], [0 0 0 0]);
HH_I = estimates(2,:);
[T2,r2] = ode45('hh', [0 499], [0 0 0 0]);
figure, hold on
plot(T1, r1(:,1), 'k')
plot(T2, r2(:,1)+120, 'k')
plot([50 150], [-40 -40], 'k')