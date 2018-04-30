% function trials = anovaTrials(bins, lambda, lambdae, ne, power) calculates
% the number of trials needed to ensure adequate power in an ANOVA to
% detect rare elevations of rate in a Poisson-process, in certain time bins. 
% 
% bins = number of time bins (equal size) per trial 
% lambda = non-elevated rate
% lambdae = elevated rate
% ne = number of elevated bins (should be << bins, e.g. 1%)
% power = proportion of experiments in which differences should be found
% (e.g. 0.8) 

function [trials, approxPower] = anovaTrials(bins, lambda, lambdae, ne, power) 
    if bins < 15 
        error('only valid for large number of bins')
    end

%     testF()

    alpha = .05;
    
    poissonVariables = 1; %use Poisson random variables (instead of Gaussian)
    
    % low precision ... 
    [trials, approxPower] = newtonsMethod(bins, lambda, lambdae, ne, power, alpha, poissonVariables, 100, [1 1000], .015); 
    
    % high precision ... 
    range = [round(trials-trials/10) round(trials+trials/10)];
    [trials, approxPower] = newtonsMethod(bins, lambda, lambdae, ne, power, alpha, poissonVariables, 1000, range, .0015);
    
    

% Newton's root-finding method 
function [trials, power] = newtonsMethod(bins, lambda, lambdae, ne, power, alpha, poissonVariables, experiments, trialRange, tolerance)

    powerRange(1) = anovaPowerExperiment(bins, lambda, lambdae, ne, trialRange(1), alpha, poissonVariables, experiments);    
    powerRange(2) = anovaPowerExperiment(bins, lambda, lambdae, ne, trialRange(2), alpha, poissonVariables, experiments);
    while powerRange(2) < power
        trialRange(2) = trialRange(2) + (trialRange(2) - trialRange(1)); %not too big a jump here if we start with a tight range
        powerRange(2) = anovaPowerExperiment(bins, lambda, lambdae, ne, trialRange(2), alpha, poissonVariables, experiments);
        trialRange
        powerRange
    end
    while powerRange(1) > power
        trialRange(1) = ceil(.75*trialRange(1)); %bigger jump OK but we don't want to go < 1
        powerRange(1) = anovaPowerExperiment(bins, lambda, lambdae, ne, trialRange(1), alpha, poissonVariables, experiments);
        trialRange
        powerRange
    end
    
    while (min(abs(powerRange - power)) > tolerance)
        t = round(mean(trialRange));
        p = anovaPowerExperiment(bins, lambda, lambdae, ne, t, alpha, poissonVariables, experiments);
        if (p - power > 0)
            trialRange = [trialRange(1) t];
            powerRange = [powerRange(1) p];
        else 
            trialRange = [t trialRange(2)];
            powerRange = [p powerRange(2)];
        end
        if (trialRange(2) - trialRange(1) <= 1)
            break;
        end
        trialRange
        powerRange
    end
    
    if abs(powerRange(1) - power) < abs(powerRange(2) - power)
        trialRange = fliplr(trialRange); 
        powerRange = fliplr(powerRange);
    end
    
    trials = trialRange(2);
    power = powerRange(2);

% % Build f distribution (tests null condition)
% function testF()
%     bins = 5;
%     trials = 20;
%     lambda = .1;
%     
%     fc = finv(.99, bins-1, (bins-1)*trials);
%     dd = fc/20;
%     domain = 0:dd:2*fc;
%     answer = fpdf(domain, bins-1, (bins-1)*trials);
%     cdf = min(cumsum(answer)*dd, 1)
%     for i = 1:length(domain)
%         approxcdf(i) = 1-anovaPowerExperiment(bins, lambda, lambda, 0, trials, 1-cdf(i), 0);      
%     end
%     approxcdf
%     approx = [0 diff(approxcdf)/dd];
%     %figure, hold on, plot(domain, answer, 'k'), plot(domain, approx, 'r')
%     figure, hold on, plot(domain, cdf, 'k'), plot(domain, approxcdf, 'r')
%     return
%         
