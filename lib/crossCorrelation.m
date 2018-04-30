% function [shift corr, peak, area, narea] = crossCorrelation(A, B, bin) calculates the
% cross-correlation between spike trains 'A' and 'B', with bin 
% size 'bin'. 
% 
% Return Args
% shift: domain fo time shifts for cross correlogram
% corr: within-bin coincidences/s at each time shift
% peak: correlation index (reported in article) from Tomita & Eggermont, 2004
% area: area of the peak region normalized by expected coincidences per bin
% narea: area of the peak region normalized by the expected area under the peak region

function [shift, corr, peak, area, narea] = crossCorrelation(A, B, binSize)
    A = A(find(A>0));
    B = B(find(B>0));    
    T = min([A(end) B(end)]);
    N = ceil(T / binSize);
    Na = length(A);
    Nb = length(B);
    E = Na * Nb / N; %total expected coincidences
    
    stepSize = .001;
    binnedA = bin(A, binSize, N);
    width = .1; %time to search forward/back for correlations (those we created are temporally localized)
    shift = -width:stepSize:width; 
    for i = 1:length(shift)
        shiftedB = B + shift(i);
        shiftedB = shiftedB(find(shiftedB > 0 & shiftedB <= T));
        binnedB = bin(shiftedB, binSize, N);
        corr(i) = sum(binnedA .* binnedB); %total coincidences at this shift
    end
    
    coeff = (corr - E) / ((Na - Na^2/N) * (Nb - Nb^2/N))^.5;
    peak = max(coeff);    

    %smoothing kernel to help find elevated area around peak
    SD = .002;
    gaussian = (1 / (SD*(2*pi)^.5)) * exp(-shift.^2/(2*SD^2));
    
    % peak area calculations
    filtered = conv(gaussian, coeff);
    strip = (length(gaussian)-1)/2;
    filtered = filtered(strip+1:end-strip)/sum(gaussian);
    peakFiltered = max(filtered);
    peakIndex = find(filtered == peakFiltered); 
    peakIndex = peakIndex(1);
    topEdge = peakIndex + findEdge(filtered(peakIndex:end), peakFiltered*.1) - 1;
    bottomEdge = peakIndex - findEdge(fliplr(filtered(1:peakIndex)), peakFiltered*.1) + 1; 
    peakArea = corr(bottomEdge:topEdge);
    area = ((sum(peakArea) - length(peakArea)*E)) / E;
    narea = ((sum(peakArea) - length(peakArea)*E)) / (length(peakArea)*E);
    
    corr = corr / T; %return as coincidences per second   
    
    if (stepSize ~= binSize) 
        warning('Step size and bin size are different: peak area result is biased (either missing some correlation or counting it multiple times)')
    end
    
% assign spikes at given times to N bins of given size
function binned = bin(spikes, binSize, N)
    binned = zeros(1, N);
    for i = 1:length(spikes)
        binNum = ceil(spikes(i) / binSize);
        if (binNum <= N)
            binned(binNum) = binned(binNum) + 1;
        end
    end

% finds the index of the edge of a correlation peak (in the forward
% direction)
function offset = findEdge(filtered, cutoff)
    indices = find(filtered < cutoff);
    if (indices)
        offset = min(indices)-1;
    else
        offset = length(filtered);
    end
    
    