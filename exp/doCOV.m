% function doCOV(mix, noise, dataFile) runs an experiment testing the
% effects of varying coefficients of variation with given parameters. 
% Called by exp_COV and exp_COVNS. 
% 
% mix: an nX3 matrix defining n COV conditions, each row containing 
%   proportions of Poisson, regular, and irregular-bursting ISI densities. 
% noise: an mX2 matrix defining m noise conditions, each row containing 
%   spike time jitter (ms) and noise spikes (proportion of non-noise
%   spikes)
% dataFile: name of file in which to save results
% useBandLimitedJitter: 0=white-noise jitter; 1=band-limited jitter

function doCOV(mix, noise, dataFile, useBandLimitedJitter)
    
    x = load('signals_figure4.mat');
    signals = x.signals;

    rate = 40;

    for i = 1:size(mix,1)
        [spikes, cov] = genUncorrelated(20, 10, .0002, rate, mix(i,:), struct('SD', .0025, 'meanSD', .0025));
        meanCOV(i) = mean(cov);
        sdCOV(i) = std(cov);

        for j = 1:size(noise,1)    

            for k = 1:size(signals,1)
                signal = signals(k,:);
                [spikes, cov] = genUncorrelated(500, .3, .0002, rate, mix(i,:), struct('SD', .0025, 'meanSD', .0025));

                jitter = noise(j,1);
                noiseRate = noise(j,2) * rate;

                nt = 32;
                ne = 5;
                if jitter == 0 & noiseRate == 0
                    nt = 1
                    ne = 1
                end
                    
                if useBandLimitedJitter
                    sprintf('Using band-limited jitter')
                    [weights, err, t] = decode(signal, .0002, spikes, [0 0 0], [jitter 0 0], [0 2*pi*15], noiseRate, nt, ne, 1); t
                else 
                    [weights, err, t] = decode(signal, .0002, spikes, [jitter 0 0], [0 0 0], [0 10], noiseRate, nt, ne, 1); t
                end
                meanErr(i,j,k) = mean(err)
                sdErr(i,j,k) = std(err);
            end
        end
        save(dataFile) 
    end

