% function w = constrainedDecoders(target, components, relNoise) calculates 
% optimal weights for approximating a target with a set of neural
% responses. The signs of the weights are constrained, so that the first
% 80% of components have positive sign, and the remainder have -ve sign. 
%
% signal: the desired composite signal
% components: must have component signals as rows 
% relNoise: noise added (proportion of max over components)

function weights = optimalPositiveWeights(target, components, relNoise)
    if (relNoise > 0)
        noise = max(max(abs(components))) * relNoise * randn(size(components));
        components = components + noise;
    end    

    global PRESENT_TARGET PRESENT_COMPONENTS
    
    PRESENT_TARGET = target;
    PRESENT_COMPONENTS = components;
    
    n = size(components,1);
    
    % lower and upper bounds on weights
    npos = floor(n * 0.8);
%     LB = [zeros(npos,1); -b*ones(n-npos,1)];
%     UB = [b*ones(npos,1); zeros(n-npos,1)];
    LB = zeros(n,1);
    UB = zeros(n,1);
    LB(npos:end) = -Inf;
    UB(1:npos) = Inf;

%     A = diag([-ones(1,npos) ones(1,n-npos)]);
%     B = zeros(n,1);
    
    % scaled inner product of each component with target as starting weights 
%     sip = components * target';
%     sip = max(0,sip);
%     sip(npos:end) = 0;
%     est = sum(components .* (sip * ones(1,length(target))),1);
%     integral = sum(max(0,target));
%     scale = integral / sum(est);
%     sip = sip * scale;
% %     figure, hold on
% %     plot(target, 'k')
% %     plot(est * scale, 'r')
% %     pause
%     x0 = sip;
    x0 = zeros(n,1);
    
    options = optimset('GradObj','on');
    
    tic
    [weights, fval] = fmincon(@weightfun,x0,[],[],[],[],LB,UB,[],options);
%     [weights, fval] = fmincon(@weightfun,x0,A,B);
    toc
    
    weighted = weights * ones(1,size(components, 2)) .* components;
    s = sum(weighted,1);
    diff = target - s;
    
    figure, hold on
    plot(target, 'k--')
    plot(s, 'k')
    plot(diff, 'r')
    pause
    
    