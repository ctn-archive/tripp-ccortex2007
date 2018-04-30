% function all = ensembleCrossCorr(spikes, binsize, maxTime) returns
% cross-correlations between all pairs of spike trains in an ensemble. 
% Correlation of each pair is given in a row, with columns: peak, area,
% narea (as defined in crossCorrelation.m)
% 
% spikes: ensemble spike times (one train per row)
% binSize: see crossCorrelation
% maxTime: time up to which trains are compared (s)

function all = ensembleCrossCorr(spikes, binSize, maxTime)
    all = [];
    for i = 1:size(spikes,1)
        for j = i+1:size(spikes,1)
            A = spikes(i,:); A = A(find(A <= maxTime & A > 0));
            B = spikes(j,:); B = B(find(B <= maxTime & B > 0));
            if (length(A) > 0 & length(B) > 0)
                [shift, corr, peak, area, narea] = crossCorrelation(A, B, binSize);
                all = [all; peak area narea];
            end
        end
    end
    