% This script experiments with subtle patterns driving Hodgkin-Huxley 
% neurons more accurately via dynamic effects

clear all

global HH_T HH_I

dt = .0002;
T = .6;
psc = PSC(dt);

% two PSC time constants will be compared 
pte = 0:dt:.2;
pte2 = 0:dt:.05;
psce = exp(-pte/.02); psce = psce/sum(psce);
psce2 = exp(-pte2/.005); psce2 = psce2/sum(psce2);

trials = 100;
for i = 1:trials
    % Assume a pool of 50,000 neurons each of which fires at a low rate over
    % the course of 1/2 second, except for one 5ms bin in which it fires at a
    % higher rate. There is an even distribution of such bins over the .5s. We
    % choose the 10,000 of these neurons with elevated bins in the 200-300 ms
    % range (again evenly distributed). Poisson firing all around. Groups of
    % 500 for each of the 20 bins from 200-300ms. 
    lowRate = 200 * .1;
    highRate = 200 * .1799; %876 trials needed with elevated rate 200*0.1799 (error: 12.2/20 spikes/s) (from exp_power results)
    n = 10000;
    groupSize = 500;
    
    % generate firing patterns for lead-in of .1s plus .5s of simulation 
    [spikes, cov] = genUncorrelated(n, T, .0002, lowRate, [1 0 0], struct('RT', 0));
    
    % generate additional elevated-rate spikes and add to above patterns at appropriate times  
    [elevSpikes, cov] = genUncorrelated(n, .05, .0002, highRate-lowRate, [1 0 0], struct('RT', 0));
    add = T - .3 + .005 * floor((0:n-1)/groupSize)' * ones(1,size(elevSpikes,2));
    elevSpikes = elevSpikes + add;
    elevSpikes(find(elevSpikes-add > .005)) = 0;
    spikes = [spikes elevSpikes]; %ordering doesn't matter here

    spikeIndices = round(spikes/dt);
    len = round(T/dt);
    spikeSignal = zeros(1,len);
    for j = 1:size(spikes,1)
        ind = spikeIndices(j,:);
        ind = ind(find(ind>0 & ind<=len));
        spikeSignal(ind) = spikeSignal(ind) + 1;
    end
    
    % two alternative excitatory inputs calculated from different PSC models 
    excitation = conv(spikeSignal, psce); excitation = excitation(1:len);
    excitation = excitation * 2;
    excitation = max(74, excitation); %clipped to reduce startup artefact

    excitation2 = conv(spikeSignal, psce2); excitation2 = excitation2(1:len);        
    excitation2 = excitation2 * 2;
    excitation2 = max(74, excitation2);

    HH_T = (dt*1000) * (1:length(excitation));
    
    % we test a range of balancing inhibition and report the best
    % results (i.e. the ones that result in the most selective firing)
    inhibition = 70:1:95;

    for j = 1:length(inhibition)
        % Hodgkin-Huxley simulations 
        HH_I = excitation - inhibition(j);        
        [odet1,r] = ode45('hh', [1 T*1000], [0 0 0 0]);
        V1 = r(:,1)';
        indices = find([0 diff(V1 > 50)] > .5);
        spikeTimes = odet1(indices) / 1000;
	
        HH_I = excitation2 - inhibition(j);
        [odet2,r] = ode45('hh', [1 T*1000], [0 0 0 0]);
        V2 = r(:,1)';
        indices = find([0 diff(V2 > 50)] > .5);
        spikeTimes2 = odet2(indices) / 1000;
        
        layer1(i,j,1:size(spikes,2)) = spikes(j,:); %example presynaptic (near-random) neuron
        layer2a(i,j,1:length(spikeTimes)) = spikeTimes; %with slow-PSC excitation 
        layer2b(i,j,1:length(spikeTimes2))= spikeTimes2; %with fast-PSC excitation 
	end
    
    save data_subtleDrive.mat    
end

