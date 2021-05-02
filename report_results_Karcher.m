%% MATLAB initials:
clc
clear all
close all

%% Settings:
experiment = "Karcher_mean";  %%--> Karcher_mean, RiemMix
retraction_type = "expm"; %--> expm , taylor
dimenion_of_matrix = 100;
n_runs = 2;
legend_of_methods = {'VTF-RLBFGS (ISR)', 'VTF-RLBFGS (Cholesky)', 'RLBFGS (Wolfe)', 'RLBFGS (Cautious)'};
iterations_to_report = {1, 2, "last"};
plot_again = true;
average_results_again = false;

%% load results
methods = {'RLBFGS_Wolfe_VTFree_'+retraction_type, 'RLBFGS_Wolfe_VTFreeCholesky_'+retraction_type, 'RLBFGS_Wolfe_'+retraction_type, 'RLBFGS_cautious'}; 
costs_all_list = cell(n_runs, length(methods));
costs_list = zeros(n_runs, length(methods), length(iterations_to_report));
n_iterations_list = zeros(n_runs, length(methods));
time_list = zeros(n_runs, length(methods));
time_average_list = zeros(n_runs, length(methods));

for method_index = 1:length(methods)
    method = methods{method_index};
    path = pwd+"/saved_files/"+experiment+"/dim="+dimenion_of_matrix+"/all_info_"+method+".mat";
    load(path);  %--> load all_info.mat
    for run_index = 1:length(all_info)
        info_ = all_info{run_index};
        all_info_list{run_index}{method_index} = info_;
        cost = convert_struct_to_array(info_, "cost");
        time = convert_struct_to_array(info_, "time");
        costs_all_list{run_index, method_index} = cost;
        for itr_of_cost_index = 1:length(iterations_to_report)
            itr_ = iterations_to_report{itr_of_cost_index};
            if isstring(itr_)
                costs_list(run_index, method_index, itr_of_cost_index) = cost(end);
            else
                if itr_ > length(cost); continue; end
                costs_list(run_index, method_index, itr_of_cost_index) = cost(itr_);
            end
            n_iterations_list(run_index, method_index) = length(cost);
            time_list(run_index, method_index) = time(end);
            time_average_list(run_index, method_index) = mean(time);
        end
    end 
end

%% plot results:
if plot_again
    for run_index = 1:length(all_info)
        PLOTFOLDER = pwd + "/saved_files/Karcher_mean/dim=" + dimenion_of_matrix + "/run" + run_index + "/plots_final/";
        if ~exist(PLOTFOLDER, 'dir')
            mkdir(PLOTFOLDER);
        end
        sim1_plot_results_forKarcherExperiment(costs_all_list, run_index, dimenion_of_matrix, PLOTFOLDER, retraction_type, legend_of_methods);
    end
end

%% take mean and std and save as table:
if average_results_again
    path_ = pwd + "/saved_files/Karcher_mean/dim=" + dimenion_of_matrix + "/";
    fid = fopen(path_+'/results_'+retraction_type+'.txt', 'wt');
    fprintf(fid, 'Method \t itr=%d \t itr=%d \t itr=last \t n_iters \t time \t average time \n', iterations_to_report{1}, iterations_to_report{2});
    fprintf(fid, '=============================================== \n');
    cost_mean_along_runs = mean(costs_list, 1); z = size(cost_mean_along_runs); cost_mean_along_runs = reshape(cost_mean_along_runs, [z(2:end) 1]);
    cost_std_along_runs = std(costs_list, 1); z = size(cost_std_along_runs); cost_std_along_runs = reshape(cost_std_along_runs, [z(2:end) 1]);
    n_itr_mean_along_runs = mean(n_iterations_list, 1); z = size(n_itr_mean_along_runs); n_itr_mean_along_runs = reshape(n_itr_mean_along_runs, [z(2:end) 1]);
    n_itr_std_along_runs = std(n_iterations_list, 1); z = size(n_itr_std_along_runs); n_itr_std_along_runs = reshape(n_itr_std_along_runs, [z(2:end) 1]);
    time_mean_along_runs = mean(time_list, 1); z = size(time_mean_along_runs); time_mean_along_runs = reshape(time_mean_along_runs, [z(2:end) 1]);
    time_std_along_runs = std(time_list, 1); z = size(time_std_along_runs); time_std_along_runs = reshape(time_std_along_runs, [z(2:end) 1]);
    time_average_mean_along_runs = mean(time_average_list, 1); z = size(time_average_mean_along_runs); time_average_mean_along_runs = reshape(time_average_mean_along_runs, [z(2:end) 1]);
    time_average_std_along_runs = std(time_average_list, 1); z = size(time_average_std_along_runs); time_average_std_along_runs = reshape(time_average_std_along_runs, [z(2:end) 1]);
    %--> save results in table format in text file:
    for method_index = 1:length(methods)
        for itr_of_cost_index = 1:length(iterations_to_report)
            mean_cost(itr_of_cost_index) = cost_mean_along_runs(method_index, itr_of_cost_index);
            std_cost(itr_of_cost_index) = cost_std_along_runs(method_index, itr_of_cost_index);
        end
        n_itr_mean = n_itr_mean_along_runs(method_index);
        n_itr_std = n_itr_std_along_runs(method_index);
        time_mean = time_mean_along_runs(method_index);
        time_std = time_std_along_runs(method_index);
        time_average_mean = time_average_mean_along_runs(method_index);
        time_average_std = time_average_std_along_runs(method_index);
        fprintf(fid, '%s \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \n', methods{method_index}, mean_cost(1), std_cost(1), mean_cost(2), std_cost(2), mean_cost(3), std_cost(3), n_itr_mean, n_itr_std, time_mean, time_std, time_average_mean, time_average_std);
    end
    fclose(fid);
end

%% functions:
function sim1_plot_results_forKarcherExperiment(costs_all_list, run_index, dim, PLOTFOLDER, retraction_type, legend_of_methods)
%---> costs_all_list --> cell(run, method)

methods = legend_of_methods;

if true
    PLOT_ALL_DIFF = true;
else
    PLOT_ALL_DIFF = false; 
end

iplot = 0;
ma = -Inf;
for imethod = 1:numel(legend_of_methods)
    cost = costs_all_list{run_index, imethod};

    iplot = iplot + 1;
    % computin max and min of y-data
    xData{iplot} = 1:length(cost);
    yData{iplot} = nanmean(cost, 1);
    if true
        ma = min(yData{iplot});
    else
        ma = max(max(yData{iplot}), ma);
    end
end

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

if false % needed for one of the methods
    pp_options.ylimits = [21.299 23.2001];
    pp_options.xlimits = [0.1 2*10^2];
end

if false
    pp_options.xlimits = [0 100];
end

yData_before_difference = yData;

if PLOT_ALL_DIFF
    if false
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
            %yData{imethod} = ma - yData{imethod};
            yData{imethod} = yData{imethod} - ma;
            %%% added by Benyamin for stability in log for plot:
            for k = 1:length(yData{imethod})
                if yData{imethod}(k) == 0
                    yData{imethod}(k) = 1e-20;
                end
            end
        end
%         %%% added by Benyamin for avoiding exact zero for stability in log for plot:
%         minimum_cost_positive = inf;
%         for imethod = 1:numel(methods)
%             for k = 1:length(yData{imethod})
%                 if yData{imethod}(k) > 0 && yData{imethod}(k) < minimum_cost_positive
%                     minimum_cost_positive = yData{imethod}(k);
%                 end
%             end
%         end
%         for imethod = 1:numel(methods)
%             for k = 1:length(yData{imethod})
%                 if yData{imethod}(k) == 0
%                     yData{imethod}(k) = minimum_cost_positive;
%                 end
%             end
%         end
%         %%%
        pp_options.ylabel = 'Averaged Cost Difference';
        pp_options.legendLoc = 'NorthEast';
    end
    pngfile = sprintf('%s/sim1_results_dim(%d)_%s.png', PLOTFOLDER, dim, retraction_type);
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
    pngfile = sprintf('%s/sim1_results_dim(%d)_%s.png', PLOTFOLDER, dim, retraction_type);
end

pp_options.logScale = LOG_SCALE;

%%%% user-defined legend added by Benyamin:
pp_options.legend = legend_of_methods;

figure('Name', pngfile, 'visible', 'off')
prettyPlot(xData, yData, pp_options);
saveas(gcf, pngfile, 'png')
%matlab2tikz(texfile);
%export_fig(pngfile, '-png')
%export_fig(pdffile, '-pdf', '-transparent')

filename = sprintf('%s/sim1_results_dim(%d)_%s.png', PLOTFOLDER, dim, retraction_type);
save(filename+"_xData.mat","xData");
save(filename+"_yData.mat","yData");
save(filename+"_yData_before_difference.mat","yData_before_difference");

end


