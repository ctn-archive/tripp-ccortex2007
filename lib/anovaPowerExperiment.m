% function power = anovaPowerExperiment(bins, lambda, lambdae, ne, trials, alpha)
% runs a numerical experiment to estimate statistical power in an ANOVA. 
% 
% bins: number of time bins per trial, among which we look for elevated firing rates
% lambda: non-elevated firing rate
% lambdae: elevated rate (assumed to occur in few bins)
% ne: number of bins that have elevated rate
% trials: number of trials in experiment
% alpha: significance level 
% 
% Optional Args 
% poissonVariables: 0=Gaussian rate variability; 1=Poisson variability (default)
% experiments: number of experiments to perform (default=1000) 
 
% Note: to check against standard tables (e.g. Cohen), use Gaussian random variables and
% change to one trial elevated and one lowered. 

function power = anovaPowerExperiment(bins, lambda, lambdae, ne, trials, alpha, varargin)
    
    poissonVariables = 1;
    experiments = 1000;
    if nargin > 6
        poissonVariables = varargin{1};
    end
    if nargin > 7
        experiments = varargin{2};
    end
    sprintf('Poisson=%i Experiments=%i', poissonVariables, experiments)
        
    cdfne = cdf('poiss', 0:20, lambda); %note: start at 0 because we use > in j-loops
    cdfe = cdf('poiss', 0:20, lambdae);
    pne = min(find(1-cdfne < .00001));
    pe = min(find(1-cdfe < .00001));
    
    fc = finv(1-alpha, bins-1, (bins-1)*trials); % critical value of ANOVA f statistic
    rejections = 0; %to keep track of # of times null hypothesis rejected 
    lambdal = lambda - (lambdae - lambda); %for Cohen table test
    for i = 1:experiments
        
        % generate random data from specified mean rates
        if (poissonVariables)
            rvalsne = rand(trials, bins-ne);
            datane = zeros(size(rvalsne));
            for j = 1:pne
                datane = datane + (rvalsne > cdfne(j));
            end        
            rvalse = rand(trials, ne);
            datae = zeros(size(rvalse));
            for j = 1:pe
                datae = datae + (rvalse > cdfe(j));
            end
            data = [datane datae];
        else
            % Cohen table test ********
            %data = [lambda + lambda^.5 * randn(trials, bins-ne)  lambdae + lambdae^.5 * randn(trials, 1)  lambdal + lambdal^.5 * randn(trials,1)]; 
            % *************************

            % Normal usage ************
            data = [lambda + lambda^.5 * randn(trials, bins-ne)  lambdae + lambdae^.5 * randn(trials, ne)];
            % *************************
        end

        % calculate ANOVA f statistic 
        xbar = mean(data);
        xbarbar = mean(xbar);
        s2 = var(data);
        sw2 = sum(s2) / bins;
        sb2 = sum( (xbar-xbarbar).^2 ) * trials / (bins-1);
        
        f = sb2/sw2;
        fs(i) = f;
        if (f > fc) 
            rejections = rejections + 1;
        end
    end
    
    power = rejections / experiments;
    
    