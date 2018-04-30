% This script generates data on error associated with signal generation by
% various Izhikevich networks

clear all

useBandLimitedJitter = 0;

cases = [ ... % distribution exponent, connectivity 
        2 .5; ...
        .1 .3; ...
        2 .65; ...
        ]
    
x = load('signals_figure4.mat');
signals = x.signals;

jitter = [0 .001 .002 .004 .008 .016];

meanPeakCorrelation = [];
meanAreaCorrelation = [];
meanNormAreaCorrelation = [];

for i = 1:size(cases,1)
    % generate long firing patterns to get statistics
    activity = izhikevichNetwork(20000, 500, cases(i,1), cases(i,2)) %note that # of neurons effects COV and synchrony (use same # as in tests)
    cov = activity.COV;
    cov = cov(find(cov > 0));
    meanCOV(i) = mean(cov);
    sdCOV(i) = std(cov);
    xc = ensembleCrossCorr(activity.spikes(251:270,:), .001, 20);
    meanPeakCorrelation = [meanPeakCorrelation; mean(xc(:,1)) std(xc(:,1))];
    meanAreaCorrelation = [meanAreaCorrelation; mean(xc(:,2)) std(xc(:,2))];
    meanNormAreaCorrelation = [meanNormAreaCorrelation; mean(xc(:,3)) std(xc(:,3))]; 

    for j = 1:length(jitter) 
        for k = 1:size(signals,1)
            % generate test firing patterns
            activity = izhikevichNetwork(300, 500, cases(i,1), cases(i,2));    
            spikes = activity.spikes;
            %figure, plot(spikes, (1:size(spikes,1))' * ones(1, size(spikes,2)), 'k.'); pause
	
            % find weights to approximate target current signal
            if useBandLimitedJitter
                sprintf('Using band-limited jitter')
                [weights, err, t] = decode(signals(k,:), .0002, spikes, [0 0 0], [jitter(j) 0 0], [0 2*pi*15], 0, 32, 5, 1); 
            else 
                [weights, err, t] = decode(signals(k,:), .0002, spikes, [jitter(j) 0 0], [0 0 0], [0 10], 0, 32, 5, 1); 
            end
            t
            meanErr(i,j,k) = mean(err)
            sdErr(i,j,k) = std(err);
        end
    end
    save 'data_Izhikevich.mat'
end

