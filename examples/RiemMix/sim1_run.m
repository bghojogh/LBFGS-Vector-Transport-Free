function sim1_run(DIM, SEP, INIT, K, N, e, RESFOLDER)

%% settings

if nargin==0
    DIM = 2; % data dimensions
    SEP = 'high';
    INIT = 'default'; % 'default', 'kmeans', 'kmeanspp'
    K = 2;
    N = 500;
end

if nargin < 6
    e = 10;
end

if nargin < 7
    RESFOLDER = 'result';
end

RUN_COUNT = 1;
% RUN_COUNT = 10;

MAXITER = 1500;
% MAXITER = 1;

DEBUG = 0;

% optimization methods
METHODS = get_methods(DIM, K, N);

% common options
co.maxiter = MAXITER;
co.tolgradnorm = -Inf;
co.minstepsize = -Inf;
co.tolcostdiff = 1e-6;
co.verbosity = 2;
co.statsfun = @statsfun;
    function stats = statsfun(D, theta, stats)
        stats.theta = theta;
        stats.ll = D.ll(theta, data) / size(data,2);
    end



%% run simulation

% load data
data_filename = sprintf('%s/sim1_data_dim(%d)_K(%d)_N(%d)_E(%d)_sep(%s)', RESFOLDER, DIM, K, N, e, SEP);
if ~exist([data_filename '.mat'], 'file')
    sim1_gen_data(DIM, K, N, e, RESFOLDER);
end
f = load(data_filename);
data = f.data;

% generate initialization params for each run
init_filename = sprintf('%s/sim1_init_dim(%d)_K(%d)_N(%d)_E(%d)_sep(%s)', RESFOLDER, DIM, K, N, e, SEP);
if ~exist([init_filename '.mat'], 'file')
    Theta0s = cell(RUN_COUNT,1);
    for run = 1:RUN_COUNT
        Theta0s{run} = init_param(data, K, INIT);
    end
    save(init_filename, 'Theta0s');
else
    fprintf('* Note: Init File exists. Skipped: %s\n', init_filename);
    load(init_filename, 'Theta0s');
end

methods = fieldnames(METHODS);
for imethod = 1:numel(methods)
    method = methods{imethod};
    
    filename = sprintf('%s/sim1_results_%s_dim(%d)_K(%d)_N(%d)_E(%d)_sep(%s)_init(%s)', RESFOLDER, method, DIM, K, N, e, SEP, INIT);
    if exist([filename '.mat'], 'file')
        fprintf('* Note: File exists. Skipped: %s\n', filename);
        continue
    end
    
    if K == 1
        D = METHODS.(method).ComponentD;
    else
        D = mixturefactory(METHODS.(method).ComponentD, K);
    end
    options = mxe_options(METHODS.(method).options, co);
    
    
    % results pre-allocation
    results = preallocate(METHODS.(method).info_fields);
    max_infoSize = 0;
    for run = 1:RUN_COUNT
        fprintf('\nMETHOD: %s, DIM: %d, SEP: %s, INIT: %s, Run: %d\n', method, DIM, SEP, INIT, run)
        
        % theta0
        theta0 = Theta0s{run};
        if K > 1
            if ~isempty(theta0)
                D2 = D.component(1);
                if isfield(D2,'main2new')
                    for k = 1:length(theta0.D)
                        theta0.D{k} = D2.main2new(theta0.D{k});
                    end
                end
            end
        else
            if isfield(D, 'main2new')
                theta0 =  D.main2new(theta0);
            end
        end
        options.theta0 = theta0;
        
        if DEBUG
            [theta, unused, info] = vizest(D, data, options, method);
        else
            [theta, unused, info] = D.estimate(data, options);
        end
        
        infoSize = numel(info);
        results(run, 1:infoSize) = info;
        max_infoSize = max(max_infoSize, infoSize);
    end
    if max_infoSize==0
        error('max_infoSize==0, %s', filename)
    end
    results = results(:, 1:max_infoSize); %#ok<NASGU>
    
    
    % save the results
    save(filename, 'results')
    
end


    function S = preallocate(fields)
        
        S(1,1) = fields;
        S(RUN_COUNT, MAXITER).iter = [];
        names = fieldnames(fields);
        for inames = 1:numel(names)
            [S.(names{inames})] = deal(nan);
        end
    end

end


