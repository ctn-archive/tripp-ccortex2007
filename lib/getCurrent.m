% function current = getCurrent(spikes, dt, len, psc) finds unscaled currents 
% arising from given spike times. Len is the the signal vector length (we 
% ignore spikes after this). 
% 
% spikes: spike times (one neuron per row)
% dt: time step of PSC and returned current
% len: # of time steps to calculate current for 
% psc: post-synaptic current kernel
% 
% current: spikes convolved with PSC (a size(spikes,1) by len matrix)

function current = getCurrent(spikes, dt, len, psc)
    spikeIndices = round(spikes/dt);
    current = zeros(size(spikes,1), len);
    for j = 1:size(spikes,1)
        ind = spikeIndices(j,:);
        ind = ind(find(ind>0 & ind<=len));
        spikeSignal = zeros(1,len);
        spikeSignal(ind) = 1;
        filtered = conv(spikeSignal, psc); %note: row-wise conv is about twice as fast as convn
        current(j,:) = filtered(:,1:len);
    end
