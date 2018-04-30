% function dy = hh(t, y) gives time derivatives of state variables for the 
% Hodgkin-Huxley model (as presented in Koch, 1999, Biophysics of
% Computation). 
%
% Uses global variables HH_T and HH_I to define external input
% (these are time and injected current -- we interpolate between given
% points). Note that injected current is in units of microamps / cm^2. 
% If starting with nanoamps, multiply by patch size (micrometres squared)
% over 100000 to get HH_I. 

% Run like this: 
% 
% global HH_T HH_I
% HH_T = 0:1:1000;
% HH_I = [... desired stimulus current (microamps / cm^2) ...]
% [T,r] = ODE45('hh', [0 1000], [0 0 0 0]);

function dy = hh(t, y)

    global HH_T HH_I
    
    % interpolate injected current
    indices = find(HH_T <= t);
    index = indices(end);
    if (HH_T(index) == t) 
        I_inj = HH_I(index);
    else 
        inTime = HH_T(index:index+1);
        inCurr = HH_I(index:index+1);
        alpha = (t - inTime(1)) / (inTime(2) - inTime(1));
        I_inj = inCurr(1) + alpha * (inCurr(2) - inCurr(1));
    end
    
    V = y(1);
    m = y(2);
    h = y(3);
    n = y(4);
    
    G_Na = 120;
    E_Na = 115; %this potential and others are relative to -60 mV
    G_K = 36;
    E_K = -12;
    G_m = 0.3;
    V_rest = 10.613;
    C_m = 1;
    
    alpha_m = (25-V) / (10 * (exp((25-V)/10) - 1));
    beta_m = 4 * exp(-V/18);
    alpha_h = 0.07 * exp(-V/20);
    beta_h = 1 / (exp((30-V)/10) + 1);
    alpha_n = (10-V) / (100 * (exp((10-V)/10) - 1));
    beta_n = 0.125 * exp(-V/80);
    
    dm = alpha_m * (1-m) - beta_m * m;
    dh = alpha_h * (1-h) - beta_h * h;
    dn = alpha_n * (1-n) - beta_n * n;
    dV = (G_Na * m^3 * h * (E_Na - V) + G_K * n^4 * (E_K - V) + G_m * (V_rest - V) + I_inj) / C_m;
    
    dy = [dV; dm; dh; dn];
    