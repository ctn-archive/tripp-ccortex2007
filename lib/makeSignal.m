% function x = makeSignal(T, dt, rms, bandwidth, randomSeed) creates a signal 
% with characteristics as specified.  
% 
% T: length of the signal (s)
% dt: time step size (s)
% rms: the desired root mean square power level of x(t)
% bandwidth: a vector with the minimum and maximum frequency components 
%   (in radians / second)
% randomSeed: a positive integer that can be used to seed randn to guarantee the generation
% of the same ("random") signal over and over again
% 
% x: resulting band-limited signal 

function x = makeSignal(T, dt, rms, bandwidth, randomSeed) 
    testPreconditions(T, dt, rms, bandwidth, randomSeed);

    if (randomSeed >= 0) 
        randn('seed', randomSeed);
    end

    omega = getFrequencies(T, dt);
    midPoint = (length(omega) + 1) / 2;

    A = zeros(size(omega));
    indices = find(omega >= bandwidth(1) & omega <= bandwidth(2));
    coeff = randn(size(indices)) + i*randn(size(indices));
    A(indices) = coeff;
    A(2*midPoint-indices) = conj(coeff);
    
    x = real(ifft(ifftshift(A)));
    [x, A] = scaleToRMS(A, x, rms);
    
% returns a list of frequency channels that can be used to compose the
% desired signal (in rad/s)
function omega = getFrequencies(T, dt) 
	dOmega = 2 * pi / T;
	n = floor(T / (2*dt));
	omega = dOmega * (-n:n);
    
% Scales the given signal so that it has the given (desired)
% root-mean-squared amplitude in time 
function [scaledSignal, scaledAmps] = scaleToRMS(amps, signal, desiredRMS)
    currentRMS = mean(signal .^ 2) .^ 0.5;
    scaleFactor = desiredRMS / currentRMS;
    scaledAmps = amps .* scaleFactor;
    scaledSignal = signal .* scaleFactor;

    
% make sure we have good parameters
function testPreconditions(T, dt, rms, bandwidth, randomSeed)
    if (dt > T) 
        error('dt should be less than T');
    end
    if (rms < 0) 
        error('rms should be zero or greater');
    end
    if (length(bandwidth) ~= 2)
        error('bandwidth should be a length-2 vector consisting of min and max allowed frequency');
    end
    if (bandwidth(1) > bandwidth(2)) 
        error('first element of bandwidth should be less than the second element');
    end
