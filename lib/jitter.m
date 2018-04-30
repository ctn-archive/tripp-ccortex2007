% function result = jitter(spikes, randSD, bandSD, band) shifts spikes 
% randomly.  
% 
% spikes: an array containing spike times (in seconds; one neuron per row)
% randSD: a 3-vector defining standard deviations for white noise (see
%   below)
% bandSD: a 3-vector defining RMS for band-limited noise
% band: a 2-vector giving low and high band limits (Hz) for band-limited
%   noise
% 
% result: as spikes but with jitter applied 
% 
% The 3-vectors above are SD for uncorrelated, correlated, and common
% jitter (in that order). Uncorrelated means that the jitter for one spike
% train is independent of jitter in the others. Common jitter is shared by all 
% neurons.  Correlated jitter is like common jitter in that it is derived
% from one source of variation for all neurons, but it is scaled by
% coefficients drawn from a distribution with mean 0 and SD 1. 

function result = jitter(spikes, randSD, bandSD, band)
    
    result = spikes;
    
    if randSD(1) > 0
        result = result + randSD(1) * randn(size(result)) .* (result > 0);
    end
    
    if max([randSD(2:3) bandSD]) > 0
        result = jitterWithDependencies(result, randSD, bandSD, band);
    end

% adds jitter that isn't independent across either time or neurons 
function result = jitterWithDependencies(spikes, randSD, bandSD, band)
    n = size(spikes,1);    
    dt = .001; %step size in which jitter is constant (< refractory time)
    time = dt:dt:max(max(spikes))+dt; 
    len = length(time);
    jitter = zeros(n, len);
    bandRMS = bandSD;
    
    corrCoeffs = randn(n,1);    
    
    % jitter functions of time applied to all neurons 
    randCorr = randSD(2) * randn(size(time));
    randComm = randSD(3) * randn(size(time));
    bandCorr = makeSignal(max(time), dt, bandRMS(2), band, -1);   bandCorr = pad(bandCorr, len);
    bandComm = makeSignal(max(time), dt, bandRMS(3), band, -1);   bandComm = pad(bandComm, len);

    % band-limited jitter unique to each neuron
    bandIndep = zeros(n, len);
    if bandRMS(1) > 0
        for i = 1:n
             bi = makeSignal(max(time), dt, bandRMS(1), band, -1); bi = pad(bi, len);
             bandIndep(i,:) = bi;
        end
    end

    % total jitter over time
    jitter =        (corrCoeffs * randCorr) + (ones(n,1) * randComm) ... 
      + bandIndep + (corrCoeffs * bandCorr) + (ones(n,1) * bandComm);
    
    
    result = zeros(size(spikes));

    % this loop is about 10% faster than the matrix method (it's the if)
    for i = 1:n
        for j = 1:size(spikes,2)
            if spikes(i,j) > 0
                result(i,j) = spikes(i,j) + jitter(i, ceil(spikes(i,j)/dt)); %last term is jitter at spike time
            end
        end
    end

% shortens or extends a band-limited signal by at most one sample to 
% reach a desired length
function result = pad(signal, desiredLength)    
    if length(signal) < desiredLength
        result = [signal signal(end)];
    elseif length(signal) > desiredLength
        result = signal(1:end-1);
    else 
        result = signal;
    end