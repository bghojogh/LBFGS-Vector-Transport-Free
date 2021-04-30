%% Settings:
experiment = "RiemMix";  %%--> Karcher_mean, RiemMix


%% load results
path = pwd+"/saved_files/"+experiment+"/all_info.mat";
load(path);



%% report results

for imethod = 1:numel(methods)
    method = methods{imethod};

    iplot = iplot+1;

    [time_run, iter_run, pll_run, ll_run, mu_mse_run, sigma_mse_run, info_] = ...
        sim1_cost_mse_likelihood_results(method, DIM, SEP, INIT, K, N, E, run, RESFOLDER);
    
    
    plot_x = iter_run;
    plot_y = pll_run;

    xData{iplot} = plot_x;
    yData{iplot} = nanmean(plot_y, 1);
end

