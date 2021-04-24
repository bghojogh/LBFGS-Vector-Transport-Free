function sim1(KIN, DIMS, SEPS, KS, NDIM, ES, INITS, SELECT)
nargin
if nargin < 1
    KIN = 11;
end

if nargin < 2
    DIMS = [2 5 20]; 
end

if nargin < 3
    % cluster separations
    SEPS = {'low', 'mid', 'high'}; %, 'no'
end

if nargin < 4
    % Number of components
    KS = [5 2]; 
end

if nargin < 5
    % N= NDIM *D^2
     NDIM = [10 100 1000]; % 100 1000
end

  
if nargin < 6
    % excentricity
    ES = [10 1];  
end

if nargin < 7
    %INITS = {'default','kmeanspp'};
    INITS = {'kmeanspp'};
end

if nargin < 8
    %INITS = {'default','kmeanspp'};
    SELECT = 'PLL';
end

addpath(genpath('/thirdparty/'))

warning off all

for KK = KIN:KIN
    RESFOLDER = sprintf('result%d',KK)
    PLOTFOLDER = sprintf('plots%d',KK)
    eval(['!mkdir ' RESFOLDER])
    eval(['!mkdir ' PLOTFOLDER])
   
    if true
        for E = ES
            for K = KS
                for DIM = DIMS
                    
                    % Number of data
                    for N = (NDIM* DIM^2)
                        for iinit = 1:numel(INITS)
                            INIT = INITS{iinit};
                            for isep = 1:numel(SEPS)
                                SEP = SEPS{isep};
                                
                                sim1_run(DIM, SEP, INIT, K, N, E, RESFOLDER)
                                
                            end
                        end
                    end
                end
                close all
            end
        end
    end
    
    if true
        for E = ES
            for K = KS
                for DIM = DIMS;
                    
                    % Number of data
                    for N = (NDIM* DIM^2)
                        for iinit = 1:numel(INITS)
                            INIT = INITS{iinit};
                            for isep = 1:numel(SEPS)
                                SEP = SEPS{isep};
                                sim1_plot_results(DIM, SEP, INIT, K, N, E, RESFOLDER, PLOTFOLDER, SELECT)
                                close all
                            end
                        end
                    end
                    
                end
            end
        end
    end
    
end
