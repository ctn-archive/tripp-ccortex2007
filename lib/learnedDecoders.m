% function [weights, errorHistory] = learnedDecoders(signal, spikes, weights, iterations, jitterSD) 
% returns learned weights with which to combine component signals in order
% to approximate a desired composite signal. 
% 
% Note: in this version, 32 distinct spike jitter results are cached and re-used
% for many iterations, which speeds things up considerably. This caching was not
% done in the published version, but it does not seem to affect the
% results. 
% 
% signal: the desired composite signal
% spikes: incoming spike pattern 
% weights: weights to start from 
% iterations: number of iterations to perform before returning 
% jitterSD: Gaussian spike jitter SD (s)
% tau: filter time constant (ignored with explicit inhibition, for
%   simplicity)
% inhibitorySpikes (optional): spike pattern for inhibitory projection (if 
%   not provided, mixed positive and negative weights are allowed as a
%   simplification)
% inhibitoryWeights (optional): initial weights for inhibitory projection

function [weights, errorHistory] = learnedDecoders(signal, spikes, weights, iterations, jitterSD, tau, varargin)

    rms = (mean(signal.^2))^.5
    
    n = size(spikes,1);
    dt = .0002;
    psc = PSC(dt);
    
    components = getCurrent(spikes, dt, length(signal), psc);

    % error signal filter  
    if tau > 0
        time = dt:dt:(tau*10);
        kernel = exp(-time/tau);
        kernel = kernel / sum(kernel);
        filtSignal = loopFilter(signal, kernel);
        filtComponents = loopFilter(components, kernel);
    end
    
    k = .00001 / max(max(components)); % learning rate ... .0001 OK for no jitter
    
    nInhibitory = 0;
    if nargin > 6
        iSpikes = varargin{1};
        nInhibitory = size(iSpikes,1);
        iComponents = getCurrent(iSpikes, dt, length(signal), psc);
        ki = .00001 / max(max(iComponents));
        if (nargin > 7) 
            iWeights = varargin{2};
        else 
            iWeights = zeros(nInhibitory,1);
        end
    end
    
    weighted = weights * ones(size(signal)) .* components;
    if nInhibitory == 0
        estimate = sum(weighted);
    else 
        iWeighted = iWeights * ones(size(signal)) .* iComponents;
        estimate = sum(weighted) + sum(iWeighted);
    end
    
    errorHistory = zeros(1,iterations+1);
    errorHistory(1) = mean( (estimate - signal).^2 );
    
    for i = 1:iterations
        if jitterSD > 0
            nt = 32;
            ne = 5;
            if i <= nt | i > iterations-ne
                jittered = jitter(spikes, [jitterSD 0 0], [0 0 0], [0 1]);
                components = getCurrent(jittered, .0002, length(signal), psc);
                if tau > 0
                    filtComponents = loopFilter(components, kernel);
                end
                if nInhibitory > 0
                    iJittered = jitter(iSpikes, [jitterSD 0 0], [0 0 0], [0 1]);
                    iComponents = getCurrent(iJittered, .0002, length(signal), psc);
                end
                if i <= nt
                    cacheComponents(:,:,i) = components;
                    if nInhibitory > 0
                        cacheIComponents(:,:,i) = iComponents;
                    end
                end
            else 
                components = cacheComponents(:,:,mod(i,nt)+1);
                if nInhibitory > 0
                    iComponents = cacheIComponents(:,:,mod(i,nt)+1);
                end
            end
        end
        
        for j = 1:length(signal)
            if nInhibitory == 0
                if tau == 0
                    E = sum(weights .* components(:,j)) - signal(j);
                    dEdw = components(:,j) * E;
                else 
                    E = sum(weights .* filtComponents(:,j)) - filtSignal(j);
                    dEdw = filtComponents(:,j) * E;
                end
                
                weights = weights - k * dEdw;
            else 
                E = sum(weights .* components(:,j)) + sum(iWeights .* iComponents(:,j)) - signal(j);
                dEdw = components(:,j) * E;
                weights = weights - k * dEdw;
                weights = max(0,weights);
                
                dEdiw = iComponents(:,j) * E;
                iWeights = iWeights - ki * dEdiw;
                iWeights = min(0, iWeights);
            end
        end
    
        weighted = weights * ones(size(signal)) .* components;
        if nInhibitory == 0
            estimate = sum(weighted);
        else 
            iWeighted = iWeights * ones(size(signal)) .* iComponents;
            estimate = sum(weighted) + sum(iWeighted);
        end
        error = mean( (estimate - signal).^2 );
        errorHistory(i+1) = error/rms;
%         sprintf('%i: %1.10f', i, error/rms)
        if mod(i,200) == 1
            i
%             figure, hold on
%             plot(signal, 'k')
%             plot(estimate, 'r');
%             pause
        end
    end

    weighted = weights * ones(size(signal)) .* components;
    if nInhibitory == 0
        estimate = sum(weighted);
    else 
        iWeighted = iWeights * ones(size(signal)) .* iComponents;
        estimate = sum(weighted) + sum(iWeighted);
    end
    
%     figure, hold on
%     plot(signal, 'k')
%     plot(estimate, 'r');
    
    if nargin > 6
        weights = [weights; iWeights];
    end
    
function result = loopFilter(signals, kernel)
    for i = 1:size(signals,1)
        filtered = conv(signals(i,:), kernel);
        trailing = filtered(size(signals,2):end);
        filtered(1:length(trailing)) = filtered(1:length(trailing)) + trailing;
        result(i,:) = filtered(1:size(signals,2));
    end
%     figure, hold on, plot(signals', 'b'), plot(result', 'k')
%     pause
    