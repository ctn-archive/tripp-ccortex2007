% function decVec = optimalDecoders(signal, components) returns 
% optimal decoding vectors with which to combine component signals in order
% to approximate a desired composite signal. 
% 
% signal: the desired composite signal
% components: must have component signals as rows 
% relNoise: noise added (proportion of max over components)

function decVec = optimalDecoders(signal, components, relNoise)
    checkDimensions(signal, components);
    
    if (relNoise > 0)
        noise = max(max(abs(components))) * relNoise * randn(size(components));
        components = components + noise;
    end
    
    gamma = components * components'; 
    invgamma = pinv(gamma);
    
    V = components * signal';
    decVec = invgamma*V;

function checkDimensions(signal, components)
    signalSize = size(signal);
    componentsSize = size(components);
    
    if (signalSize(1) ~= 1) 
        error('only 1 dimension is supported here')
    end
    
    if (signalSize(2) ~= componentsSize(2))
        error(sprintf('components length (%i) and signal length (%i) must be the same', componentsSize(2), signalSize(2)))
    end
