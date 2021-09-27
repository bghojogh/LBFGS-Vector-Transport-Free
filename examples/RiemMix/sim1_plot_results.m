function infos = sim1_plot_results(DIM, SEP, INIT, K, N, E, RESFOLDER, PLOTFOLDER, SELECT, legend_of_methods)


if nargin < 5
    N = DIM^2*100;
end

if nargin < 6
    E = 10;
end

if nargin < 7
    RESFOLDER = 'result';
end

if nargin < 8
    PLOTFOLDER = 'plots';
end

if nargin < 9
    SELECT = 'PLL';
end

PLOT_BY_TIME = false;
PLOT_TIME_DIVISIONS = 100;

if strcmp(SELECT, 'PLL')
    PLOT_ALL_DIFF = true;
else
    PLOT_ALL_DIFF = false; 
end
  
PLOTNAME = SELECT;



% optimization methods
METHODS = get_methods(DIM, K);

% run to plot
run = 1; % this can go to d1 

iplot = 0;
ma = -Inf;
methods = fieldnames(METHODS);
for imethod = 1:numel(methods)
    method = methods{imethod};
    
    [time_run, iter_run, pll_run, ll_run, mu_mse_run, sigma_mse_run, info_] = ...
        sim1_cost_mse_likelihood_results(method, DIM, SEP, INIT, K, N, E, run, RESFOLDER);
    infos(imethod).info=info_;
    infos(imethod).method = method;
    
    iplot = iplot+1;
    
    if PLOT_BY_TIME
        % interpolate results on equally-spaced time divisions
        %max_time = max(time_run(run,d2));
        max_time = max(time_run(run,:));
        plot_x = linspace(0, max_time, PLOT_TIME_DIVISIONS);
        %plot_y = nan(d1, PLOT_TIME_DIVISIONS);
        plot_y = nan(1, PLOT_TIME_DIVISIONS);
        idx_inside = (plot_x <= time_run(end));
        plot_x_run = plot_x(idx_inside);
        plot_y(run,1:numel(plot_x_run)) = interp1(time_run, ...
            ll_run, plot_x_run, 'linear'); %spline
        plot_y(run,1) = ll_run(1);
        pp_options.xlabel = 'Time (seconds)';
    else
        plot_x = iter_run;
        plot_y = pll_run;
        pp_options.xlabel = 'Iterations';
    end
    % computin max and min of y-data
    if strcmp(SELECT, 'MU')
        plot_y = mu_mse_run;
    elseif strcmp(SELECT, 'SIGMA')
        plot_y = sigma_mse_run;
    end
    xData{iplot} = plot_x;
    yData{iplot} = nanmean(plot_y, 1);
    if ~strcmp(SELECT, 'PLL')
        ma = min(yData{iplot});
    else
        ma = max(max(yData{iplot}), ma);
    end
    pp_options.legend{1,iplot} = METHODS.(method).legend;
end


%
if false
    LOG_SCALE = 0;
    pp_options.colors = [0.5 0 0; 0 0 .5; 0 .5 0; .5 0 0; 0 .5 0; 0 .5 0];
    pp_options.lineStyles = {'-', '-', '-'}; %-- % :
else
    % color-style comparing just reparametrized version and usual
    LOG_SCALE = 1;
    pp_options.colors = [0.5 0 0; 0 0 .5; 0.5 .3 0; 0 0.3 0.5; 0.5 0.5 0; ...
        0.5 0 0.5; 0.3 0.5 0;0 0.5 0.3;0.5 0 0.3;0.3 0 0.5; ...
        0.5 0.15 0.15; 0.15 0.5 0.15; 0.15 0.15 0.5]*2;
    pp_options.colors = distinguishable_colors(numel(methods));
end

%
if false % needed for one of the methods
    pp_options.ylimits = [21.299 23.2001];
    pp_options.xlimits = [0.1 2*10^2];
end

if false
    pp_options.xlimits = [0 100];
end

yData_before_difference = yData;


%
if PLOT_ALL_DIFF
    if ~strcmp(SELECT, 'PLL')
        for imethod = 1:numel(methods)
            yData{imethod} = yData{imethod} - ma;
        end
        if strcmp(SELECT, 'MU')
            pp_options.ylabel = 'MSE Mean Difference';
        else
            pp_options.ylabel = 'MSE Covariance Difference';
        end
        pp_options.legendLoc = 'NorthEast';
    else
        for imethod = 1:numel(methods)
            yData{imethod} = ma - yData{imethod};
        end
        %pp_options.ylabel = 'Averaged Cost Difference';
        pp_options.ylabel = 'Cost Difference';
        pp_options.legendLoc = 'NorthEast';
    end
    % pngfile = sprintf('%s/sim1_diff_dim(%d)_K(%d)_sep(%s)_init(%s)_%s.png', PLOTFOLDER, DIM, K, SEP, INIT, PLOTNAME);
    pngfile = sprintf('%s/sim1_results_dim(%d)_K(%d)_N(%d)_E(%d)_sep(%s)_init(%s)_%s.png', PLOTFOLDER, DIM, K, N, E, SEP, INIT, PLOTNAME);  %--> Benyamin added a more complete name so plots do not overwrite each other
    LOG_SCALE = 2;
    pp_options.ylimits = [10^-5 inf];
else
    LOG_SCALE = 0;
    if ~strcmp(SELECT, 'PLL')
        pp_options.legendLoc = 'NorthEast';
        if strcmp(SELECT, 'MU')
            pp_options.ylabel = 'MSE of Mean';
        else
            pp_options.ylabel = 'MSE of Covariance';
        end
    else
        pp_options.legendLoc = 'SouthEast';
        pp_options.ylabel = 'Averaged Cost';
    end
    pngfile = sprintf('%s/sim1_results_dim(%d)_K(%d)_N(%d)_E(%d)_sep(%s)_init(%s)_%s.png', PLOTFOLDER, DIM, K, N, E, SEP, INIT, PLOTNAME);
end

pp_options.logScale = LOG_SCALE;

%%%% user-defined legend added by Benyamin:
if nargin == 10
    pp_options.legend = legend_of_methods;
end

figure('Name', pngfile, 'visible', 'off')
prettyPlot(xData, yData, pp_options);
saveas(gcf, pngfile, 'png')
%matlab2tikz(texfile);
%export_fig(pngfile, '-png')
%export_fig(pdffile, '-pdf', '-transparent')

filename = sprintf('%s/sim1_results_dim(%d)_K(%d)_N(%d)_E(%d)_sep(%s)_init(%s)_%s', PLOTFOLDER, DIM, K, N, E, SEP, INIT, PLOTNAME);
save(filename+"_xData.mat","xData");
save(filename+"_yData.mat","yData");
save(filename+"_yData_before_difference.mat","yData_before_difference");

