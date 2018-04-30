% function error = frequency(spikes, frequencies, jitter) reports error in
% decoding sine waves at the given frequencies. 
% 
% spikes: a set of spike trains (in rows) to decode 
% frequencies: frequencies to test (Hz)
% jitter: spike jitter SD (s) with which to corrupt spike trains
% 
% error: matrix of mean-squared errors with sinusoids of RMS=1 
%   (each row corresponds to a given frequency, each column to a different 
%   test phase at that frequency)

function error = frequency(spikes, frequencies, jitter) 

    dt = .0002;
    T = .5;
    time = dt:dt:T;
    
    if jitter == 0
        nt = 1;
        ne = 1;
        interval = 1;
    else 
        nt = 32;
        ne = 5;
        interval = 2;
    end
    
    for i = 1:length(frequencies)
        phase = (.2:.2:1)*2*pi;
        for j = 1:length(phase)
            tic
            signal = (1/.7071) * sin(2*pi*frequencies(i)*time + phase(j)); % RMS=1
            [weights, err, t] = decode(signal, dt, spikes, [jitter 0 0], [0 0 0], [0 10], 0, nt, ne, interval, 0);
            error(i,j) = mean(err);
            toc
        end
    end
    