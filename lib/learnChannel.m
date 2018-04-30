% function outSpikes = learnChannel(inSpikes, n, convergence) learns a
% communication channel between an ensemble with a given spike pattern and 
% a new output ensemble of LIF neurons, and returns the spikes of the 
% output ensemble. 
% 
% inSpikes: spikes of the presynaptic ensemble
% n: size of the post-synaptic LIF ensemble
% convergence: number of pre-synaptic neurons connected to each
%   post-synaptic neuron
% iterations: training iterations over data 

function outSpikes = learnChannel(inSpikes, n, convergence, iterations) 

    psc = PSC(.0002);
    current = getCurrent(inSpikes, .0002, 1500, psc);
    meanCurrent = mean(current,2);

    rate = .00000000005; % learning rate
    
    tauRef = .002;
    tauRC = .02;
    intercept = 0 + rand(n, 1) * .9
    maxFR = 200 + rand(n,1) * 200;
	x = 1 ./ (1 - exp( (tauRef - (1 ./ maxFR)) / tauRC));
	scale = (x - 1) ./ (1 - intercept)	
	bias = 1 - scale .* intercept

    alifInput = [];
    for i = 1:n
        tic
        indices = randomIndices(size(inSpikes,1), convergence);
        magnitude = mean(sum(current(indices,:)));    
        weights = rand(length(indices),1) * 1/magnitude;        
        runningMean = 0;
        maxWeight = 1/magnitude;
        
        weighted = weights * ones(1,size(current,2)) .* current(indices,:);
%         figure, hold on
%         plot(sum(weighted), 'b');
        
        weightHistory = zeros(length(indices), iterations*size(current,2));
        
        for j = 1:iterations            
%             tic
            for k = 1:size(current, 2)
                weighted = weights .* current(indices,k);
                in = bias(i) + scale(i)*sum(weighted);
                
                if (in > 1)
                    b = 1 / ( tauRef - tauRC * log(1 - 1/in) );
                else 
                    b = 0;
                end

                nMean = (j-1)*size(current,2) + k;
                runningMean = ( (nMean-1)*runningMean + b ) / nMean;
                
%                 dw = rate * (current(indices,k) - meanCurrent(indices)) * (b - runningMean) * (b > 0) * (b < maxFR(i));
                dw = rate * (current(indices,k) - meanCurrent(indices)) * (b - runningMean);
%                 sprintf('%2.15f, %2.15f, %f, %f', min(dw), max(dw), b, b - runningMean)
                weights = weights + dw;
                weights = min(maxWeight, max(0, weights));

                weightHistory(:,nMean) = weights;
                
            end
%             toc
        end
        
        weighted = weights * ones(1,size(current,2)) .* current(indices,:);
        alifInput = [alifInput; sum(weighted)];
%         plot(sum(weighted), 'r.');
%         figure, plot(weightHistory')
       toc 
    end
    
    [outSpikes, firings] = ALIF(.0002, alifInput, bias, scale, ones(size(bias))*tauRC, ones(size(bias))*tauRef, ones(size(bias)), zeros(size(bias)));


    
function indices = randomIndices(n, convergence)
    indices = 1:n;
    while length(indices) > convergence
        remove = ceil(rand * length(indices));
        indices = [indices(1:remove-1) indices(remove+1:end)];
    end
    
        
    