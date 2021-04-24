addpath(genpath('./thirdparty/'))
DIMS = [2]; % Dimension
SEPS = {'low','mid','high'}; % Separation
KS = [2]; % Number of Components
%NDIM = [10 100 1000]; % Number of Data = NDIM*DIM^2
NDIM = [10]; % Number of Data = NDIM*DIM^2
ES = [10]; % Eccentricity
INITS = {'kmeanspp'}; % Initialization
iter_End = 2; % Number of Runs
% SELECT
SELECT = 'PLL'; % 'SIGMA', 'PLL', 'MU'
%Run 2 Different Runs of Algorithm
for KIN = 1:iter_End
    sim1(KIN, DIMS, SEPS, KS, NDIM, ES, INITS, SELECT);
end
sim1_gen_table(iter_End, DIMS(1), SEPS, KS(1), NDIM, ES(1), INITS{1}, SELECT);
