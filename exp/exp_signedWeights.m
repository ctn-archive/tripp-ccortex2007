% This script experiments with sensitivity to spike jitter with varying ISI
% coefficients of variation, as exp_COV, but with weights of predetermined 
% sign.   

useBandLimitedJitter = 0;

mix = [0 1 0; 1 8 0; 1 1 0; 1 0 0; 3 0 1; 1 0 1; 1 0 5]; %weights for combining Poisson, constant, and bursting ISI PDFs 
% mix = [1 0 0; 3 0 1]; %weights for combining Poisson, constant, and bursting ISI PDFs 
noise = [0 0; .004 0];
% noise = [.004 0];

doCOV(mix, noise, 'data_signedWeights3', useBandLimitedJitter)
