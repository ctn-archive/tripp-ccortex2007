% function [spikes, firingRate] = genCoarseCorrelated(n, dt, signal, threshold, rate)
% produces spikes trains with rates that vary with the instantaneous value of a 
% template signal.  
% 
% n: number of spike trains
% dt: time step of template signal and spike generation 
% signal: the template signal
% threshold: value of template signal that corresponds to zero firing rate
% rate: mean firing rate over course of signal
% 
% firingRate: firing rate over time, derived from signal and threshold args

function [spikes, firingRate] = genCoarseCorrelated(n, dt, signal, threshold, rate)
    time = dt:dt:dt*length(signal);
    firingContour = max(0, (signal-threshold));
    firingRate = (firingContour / sum(firingContour)) * rate * time(end) / dt;
    
    spikes = zeros(n, rate*2); 
    ISI = zeros(n, rate*2); 
    for i = 1:n
        spikesi = doCoarseCorrelated(dt, firingRate, time, signal);
        spikes(i,1:length(spikesi)) = spikesi;
        
        isi = diff(spikesi);
        ISI(i,1:length(isi)) = isi;
    end

function spikes = doCoarseCorrelated(dt, firingRate, time, signal)
    spikes = [];
    for i = 1:length(time)
        instantRate = firingRate(i);
        probSpike = instantRate * dt;
        if (probSpike > rand) 
            spikes = [spikes, time(i)];
        end
    end
        
