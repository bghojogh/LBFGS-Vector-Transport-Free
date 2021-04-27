function [time_run, iter_run, pll_run, ll_run, mu_mse_run, sigma_mse_run, info_] = ...
    sim1_cost_mse_likelihood_results(method, DIM, SEP, INIT, K, N, E, run, RESFOLDER)
%
%
% This function calculates likelihood, penalized cost, MSEs, times and iterations from the saved results
%

% Specify which run from initializatino should be returned
% Sice one initialization is used d1 is equal to one
if nargin < 8
    run = 1;
end

if nargin < 9
    RESFOLDER = 'result';
end

METHODS = get_methods(DIM, K);

filename = sprintf('%s/sim1_results_%s_dim(%d)_K(%d)_N(%d)_E(%d)_sep(%s)_init(%s)', RESFOLDER, method, DIM, K, N, E, SEP, INIT);
f = load(filename);
[d1,d2] = size([f.results.cost]);

if isfield (f.results(1, 1).theta.D{1},'sigmat') && true
    for r = 1:d1
        results = f.results;
        D = mvn2factory(size(results(r,1).theta.D{1}.sigmat,1)-1);
        for k = 1:length(results)
            for k2 = 1:length(results(r,1).theta.D)
                if isfield(results(r,k).theta, 'D')
                    results(r,k).theta.D{k2} = D.main2new(D.new2main(results(r,k).theta.D{k2}));
                end
            end
        end
        D2 = mixturefactory(D, length(results(r,1).theta.D));
        data_filename = sprintf('%s/sim1_data_dim(%d)_K(%d)_N(%d)_E(%d)_sep(%s)', RESFOLDER, DIM, K, N, E, SEP);
        load(data_filename)
        for k = 1:length(results)
            %if ~isfield(METHODS.(method).options, 'penalize')
            %    error('No Penalization');
            %end
            options = mxe_options(METHODS.(method).options);
            options.penalizertheta = D2.penalizerparam(data);
            if isfield( results(k).theta, 'D' )
                cost = mxe_costgrad( D2, results(r,k).theta, data, options );
                f.results(r,k).cost = cost;
            end
        end
        results = f.results;
        %save(filename, 'results');
    end
end

% tiny problem in saved files
f.results = f.results(:,1:end-1);

% Benyamin added for saving info:
info_ = f.results;

[d1,d2] = size([f.results.cost]);
pll_all = reshape(-[f.results.cost], [d1,d2]);
time_all = reshape([f.results.time], [d1,d2]);
ll_all = reshape([f.results.ll], [d1,d2]);


% Calculating MSE distance between true and estimated parameters
data_filename = sprintf('%s/sim1_data_dim(%d)_K(%d)_N(%d)_E(%d)_sep(%s)', RESFOLDER, DIM, K, N, E, SEP);
load(data_filename, 'theta_gen');
mu_mse_all = zeros(d1, d2);
sigma_mse_all = zeros(d1, d2);
for r = 1:d1
    for k = 1:d2
        if isfield(f.results(r, k).theta, 'D' )
            [mu_dists, sigma_dists] = mseDist(f.results(r, k).theta, theta_gen);
            [C, mu_mse_all(r, k)] = hungarian(mu_dists); %.^2
            [C, sigma_mse_all(r, k)] = hungarian(sigma_dists); %.^2
            %mu_mse_all(r, k)
            %sigma_mse_all(r, k)
            %pasue
        end
    end
end

time_run = time_all(run,:);
ll_run = ll_all(run,:);
pll_run = pll_all(run,:);
mu_mse_run = mu_mse_all(run,:);
sigma_mse_run = sigma_mse_all(run,:);
nans = isnan(time_run);
time_run(nans) = [];
ll_run(nans) = [];
pll_run(nans) = [];
mu_mse_run(nans) = [];
sigma_mse_run(nans) = [];
d2_run = length(pll_run);
if numel(ll_run) ==1
    ll_run(2) = ll_run + eps;
    time_run(2) = time_run + eps;
end
iter_run = 0:d2_run-1; % For EM Algorithm
if strncmpi(method,'SGD',3)
    iter_run = 0:d2_run-1;
end
if isfield(f.results(1),'linesearch')
    iter_run(1)=0;
    for k=2:d2
        if strncmpi(method,'CG',2)
            % since only one gradient computation is needed
            % for adaptive linesearch
            evals = (f.results(run, k).linesearch.costevals + 1)/2;
        else
            evals = f.results(run, k).linesearch.costevals;
        end
        iter_run(k)= iter_run(k-1)+ evals;
    end
end
