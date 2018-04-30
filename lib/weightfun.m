% function f = weightfun(x) evaluates the objective function for 
% constrained decoder optimization. x is a vector of decoding weights. 
% 
% Relies on global variables: 
% PRESENT_TARGET: target function/vector (1 x length of function)
% PRESENT_COMPONENTS: matrix with neuron response functions as rows (N x
% length of function)

function [f, g] = weightfun(x)
    
    global PRESENT_TARGET PRESENT_COMPONENTS I
    
    I = I+1;
    
    weighted = x * ones(1,size(PRESENT_COMPONENTS, 2)) .* PRESENT_COMPONENTS;
    diff = PRESENT_TARGET - sum(weighted,1);
    
    if mod(I,1000) == 0 
        figure, hold on, plot(PRESENT_TARGET, 'k'), plot(sum(weighted,1), 'r')
    end
    
    f = 0.5 * sum(diff.^2);
    mean(diff.^2)
    
    if nargout > 1
        g = -.1 * PRESENT_COMPONENTS * diff';
    end
    