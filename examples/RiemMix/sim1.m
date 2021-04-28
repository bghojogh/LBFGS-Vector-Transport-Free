function info_list = sim1(KIN, DIMS, SEPS, KS, NDIM, ES, INITS, SELECT, path_save)
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
%     RESFOLDER = sprintf('result%d',KK);
%     PLOTFOLDER = sprintf('plots%d',KK);
%     eval(['!mkdir ' RESFOLDER])
%     eval(['!mkdir ' PLOTFOLDER])
    RESFOLDER = sprintf("%sresult",path_save);
    PLOTFOLDER = sprintf("%splots",path_save);
    if ~exist(RESFOLDER, 'dir')
        mkdir(RESFOLDER);
    end
    if ~exist(PLOTFOLDER, 'dir')
        mkdir(PLOTFOLDER);
    end
   
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
    
    counter = 1;
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
                                info_list(counter).info_list = sim1_plot_results(DIM, SEP, INIT, K, N, E, RESFOLDER, PLOTFOLDER, SELECT);
                                info_list_(counter).eccentricity = E;
                                info_list(counter).components = K;
                                info_list(counter).dim = DIM;
                                info_list(counter).separation = SEP;
                                info_list(counter).n_data = N;
                                info_list(counter).select_mode=SELECT;                                
                                info_list(counter).iinit=INIT;                                
                                info_list(counter).result_folder=RESFOLDER;                                
                                info_list(counter).plot_folder=PLOTFOLDER;                                
                                counter = counter + 1;
                                close all
                            end
                        end
                    end
                    
                end
            end
        end
    end
    
end
