%% MATLAB initials:
clc
clear all
close all

%% installing for adding paths:
addpath(genpath(fullfile("./", 'utils')))
install;

%% Settings:
experiment = "RiemMix";  %%--> Karcher_mean, RiemMix --> it should be RiemMix in this file
plot_again = false;
average_results_again = true;

%% load results
path = pwd+"/saved_files/"+experiment+"/all_info.mat";
load(path);  %--> load all_info.mat

%% plot results
%%%%%%% Note: You can set some of the methods "false" in the function get_methods.m (in path ./examples/RiemMax/) to exclude some of the methods in plots
if plot_again
    %legend_of_methods = {'VTF-RLBFGS (ISR)', 'VTF-RLBFGS (Cholesky)', 'RLBFGS (Wolfe)', 'RLBFGS (Cautious)'};
    legend_of_methods = {'VTF-RLBFGS (ISR)', 'VTF-RLBFGS (Cholesky)', 'RLBFGS', 'RLBFGS (Cautious)'};
    for experiment_index = 1:length(all_info)
        info_ = all_info(experiment_index);
        DIM = info_.dim;
        SEP = info_.separation; 
        INIT = info_.iinit;
        K = info_.components;
        N = info_.n_data;
        E = [10]; % Eccentricity
        SELECT = info_.select_mode;
        RESFOLDER = info_.result_folder;
        %PLOTFOLDER = info_.plot_folder;
        tmp = info_.plot_folder.split("plots");
        PLOTFOLDER = tmp(1) + "plots_final2";  %--> tmp(1) + "plots_final"; or tmp(1) + "plots_final2"; or ... --> can rename not to overwrite previous plots
        if ~exist(PLOTFOLDER, 'dir')
           mkdir(PLOTFOLDER)
        end
        sim1_plot_results(DIM, SEP, INIT, K, N, E, RESFOLDER, PLOTFOLDER, SELECT, legend_of_methods);
    end
end

%% average results
SEPS = {'low','mid','high'}; % Separation
methods = {'LBFGS1', 'LBFGS2', 'LBFGS3', 'LBFGS4'};
N_list = {40};  %--> {40, 400}; --> 100, 1000, 10000, 1000000 ---> it depends on the dimensionality of data --> it is 100*(dim^2)
n_runs = 10;
iterations_to_report = {30, 50, "last"};  %--> {10, 20, "last"}, {20, 50, "last"}, {30, 50, "last"}
if average_results_again
    costs_list = zeros(n_runs, length(SEPS), length(methods), length(N_list), length(iterations_to_report));
    n_iterations_list = zeros(n_runs, length(SEPS), length(methods), length(N_list));
    time_list = zeros(n_runs, length(SEPS), length(methods), length(N_list));
    time_average_list = zeros(n_runs, length(SEPS), length(methods), length(N_list));
    for experiment_index = 1:length(all_info)
        info_ = all_info(experiment_index);
%         for run_index = 1:n_runs
            aa = info_.run.split("run");
            RUN = str2double(aa(2));
            for SEP_index = 1:length(SEPS)
                SEP = SEPS{SEP_index};
                for method_index = 1:length(methods)
                    method = methods{method_index};
                    for n_data_index = 1:length(N_list)
                        N = N_list{n_data_index};
                        %%%% information:
                        DIM = info_.dim;
                        %SEP = info_.separation; 
                        INIT = info_.iinit;
                        K = info_.components;
                        %N = info_.n_data;
                        E = [10]; % Eccentricity
                        SELECT = info_.select_mode;
                        %%%% results:
                        RESFOLDER = info_.result_folder;
                        tmp = RESFOLDER.split("result");
                        RESFOLDER_final = tmp(1) + "result_final";
                        if ~exist(RESFOLDER_final, 'dir')
                           mkdir(RESFOLDER_final)
                        end
                        %%%% load result:
                        filename = sprintf('%s/sim1_results_%s_dim(%d)_K(%d)_N(%d)_E(%d)_sep(%s)_init(%s).mat', RESFOLDER, method, DIM, K, N, E, SEP, INIT);
                        load(filename);
                        cost = convert_struct_to_array(results, "cost");
                        time = convert_struct_to_array(results, "time");
                        for itr_of_cost_index = 1:length(iterations_to_report)
                            itr_ = iterations_to_report{itr_of_cost_index};
                            if isstring(itr_)
                                costs_list(RUN, SEP_index, method_index, n_data_index, itr_of_cost_index) = cost(end);
                            else
                                if itr_ > length(cost); continue; end
                                costs_list(RUN, SEP_index, method_index, n_data_index, itr_of_cost_index) = cost(itr_);
                            end
                        end
                        n_iterations_list(RUN, SEP_index, method_index, n_data_index) = length(cost);
                        time_list(RUN, SEP_index, method_index, n_data_index) = sum(time);
                        time_average_list(RUN, SEP_index, method_index, n_data_index) = mean(time);
                    end
                end
            end
%         end
    end
    %--> take mean and std:
    path_ = RESFOLDER.split("run");
    fid = fopen(path_(1)+'/results.txt', 'wt');
    fprintf(fid, 'SEP \t Method \t N \t itr=%d \t itr=%d \t itr=last \t n_iters \t time \t average time \n', iterations_to_report{1}, iterations_to_report{2});
    fprintf(fid, '=============================================== \n');
    %--> save results in table format in text file:
    for SEP_index = 1:length(SEPS)
        for method_index = 1:length(methods)
            for n_data_index = 1:length(N_list)
                for itr_of_cost_index = 1:length(iterations_to_report)
                    mean_cost(itr_of_cost_index) = mean(costs_list(:, SEP_index, method_index, n_data_index, itr_of_cost_index));
                    std_cost(itr_of_cost_index) = std(costs_list(:, SEP_index, method_index, n_data_index, itr_of_cost_index));
                end
                n_itr_mean = mean(n_iterations_list(:, SEP_index, method_index, n_data_index));
                n_itr_std = std(n_iterations_list(:, SEP_index, method_index, n_data_index));
                time_mean = mean(time_list(:, SEP_index, method_index, n_data_index));
                time_std = std(time_list(:, SEP_index, method_index, n_data_index));
                time_average_mean = mean(time_average_list(:, SEP_index, method_index, n_data_index));
                time_average_std = std(time_average_list(:, SEP_index, method_index, n_data_index));
                fprintf(fid, '%s \t %s \t %d \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \n', SEPS{SEP_index}, methods{method_index}, N_list{n_data_index}, mean_cost(1), std_cost(1), mean_cost(2), std_cost(2), mean_cost(3), std_cost(3), n_itr_mean, n_itr_std, time_mean, time_std, time_average_mean, time_average_std);
            end
        end
    end
    fclose(fid);
end




