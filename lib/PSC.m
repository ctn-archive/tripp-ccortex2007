% function PSC = PSC(dt) returns a stereotyped post-synaptic current 
% kernel (adapted from Destexhe et al, 1998, in Koch & Segev, Methods in 
% Neuronal Modelling).  

function PSC = PSC(dt)
    a = 1100000; %M^-1 s^-1
    B = 190; %s^-1
    T = .001; %1mM (as figure 4 in Destexhe et al)
    gAMPA = 1; %this is irrelevant, as we normalize below and scale by synaptic weight separately 
    V = -65;  %mV, assuming constant, far from reversal potential (-65 is a typical reversal for leak current)
    EAMPA = 0; %mV
    
    pulseTime = dt:dt:.001;
    decayTime = dt:dt:.099;
    
    ra = -a*T/(a*T+B) * exp(-(a*T+B)*pulseTime) + a*T/(a*T+B);
    rb = ra(end) * exp(-(B)*decayTime);
    
    PSC = [ra rb] * gAMPA * (EAMPA - V); %interpreting inward current as +ve
    PSC = PSC / (sum(PSC)*dt);