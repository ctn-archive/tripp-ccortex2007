% This script experiments with sensitivity to spike jitter with varying ISI
% coefficients of variation 

useBandLimitedJitter = 0;

mix = [0 1 0; 1 8 0; 1 1 0; 1 0 0; 3 0 1; 1 0 1; 1 0 5]; %weights for combining Poisson, constant, and bursting ISI PDFs 
noise = [0 0; .001 0; .002 0; .004 0; .008 0; .016 0];

doCOV(mix, noise, 'data_COV', useBandLimitedJitter)
