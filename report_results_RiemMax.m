%% MATLAB initials:
clc
clear all
close all

%% Settings:
experiment = "RiemMix";  %%--> Karcher_mean, RiemMix
plot_again = false;
average_results_again = true;

%% load results
path = pwd+"/saved_files/"+experiment+"/all_info.mat";
load(path);  %--> load all_info.mat

%% plot results
if plot_again
    legend_of_methods = {'VTF-RLBFGS (ISR)', 'VTF-RLBFGS (Cholesky)', 'RLBFGS (Wolfe)', 'RLBFGS (Cautious)'};
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
        PLOTFOLDER = tmp(1) + "plots_final";
        if ~exist(PLOTFOLDER, 'dir')
           mkdir(PLOTFOLDER)
        end
        sim1_plot_results(DIM, SEP, INIT, K, N, E, RESFOLDER, PLOTFOLDER, SELECT, legend_of_methods);
    end
end

%% average results
SEPS = {'low','mid','high'}; % Separation
methods = {'LBFGS1', 'LBFGS2', 'LBFGS3', 'LBFGS4'};
N_list = {40, 400};
n_runs = 2;
iterations_to_report = {10, 25, "last"};
if average_results_again
    costs_list = zeros(n_runs, length(SEPS), length(methods), length(N_list), length(iterations_to_report));
    n_iterations_list = zeros(n_runs, length(SEPS), length(methods), length(N_list));
    time_list = zeros(n_runs, length(SEPS), length(methods), length(N_list));
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
                        time_list(RUN, SEP_index, method_index, n_data_index) = time(end);
                    end
                end
            end
%         end
    end
    %--> take mean and std:
    path_ = RESFOLDER.split("run");
    fid = fopen(path_(1)+'/results.txt', 'wt');
    fprintf(fid, 'SEP \t Method \t N \t itr=%d \t itr=%d \t itr=last \t n_iters \t time \n', iterations_to_report{1}, iterations_to_report{2});
    fprintf(fid, '=============================================== \n');
    cost_mean_along_runs = mean(costs_list, 1); z = size(cost_mean_along_runs); cost_mean_along_runs = reshape(cost_mean_along_runs, [z(2:end) 1]);
    cost_std_along_runs = std(costs_list, 1); z = size(cost_std_along_runs); cost_std_along_runs = reshape(cost_std_along_runs, [z(2:end) 1]);
    n_itr_mean_along_runs = mean(n_iterations_list, 1); z = size(n_itr_mean_along_runs); n_itr_mean_along_runs = reshape(n_itr_mean_along_runs, [z(2:end) 1]);
    n_itr_std_along_runs = std(n_iterations_list, 1); z = size(n_itr_std_along_runs); n_itr_std_along_runs = reshape(n_itr_std_along_runs, [z(2:end) 1]);
    time_mean_along_runs = mean(time_list, 1); z = size(time_mean_along_runs); time_mean_along_runs = reshape(time_mean_along_runs, [z(2:end) 1]);
    time_std_along_runs = std(time_list, 1); z = size(time_std_along_runs); time_std_along_runs = reshape(time_std_along_runs, [z(2:end) 1]);
    %--> save results in table format in text file:
    for SEP_index = 1:length(SEPS)
        for method_index = 1:length(methods)
            for n_data_index = 1:length(N_list)
                for itr_of_cost_index = 1:length(iterations_to_report)
                    mean_cost(itr_of_cost_index) = cost_mean_along_runs(SEP_index, method_index, n_data_index, itr_of_cost_index);
                    std_cost(itr_of_cost_index) = cost_std_along_runs(SEP_index, method_index, n_data_index, itr_of_cost_index);
                end
                n_itr_mean = n_itr_mean_along_runs(SEP_index, method_index, n_data_index);
                n_itr_std = n_itr_std_along_runs(SEP_index, method_index, n_data_index);
                time_mean = time_mean_along_runs(SEP_index, method_index, n_data_index);
                time_std = time_std_along_runs(SEP_index, method_index, n_data_index);
                fprintf(fid, '%s \t %s \t %d \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \n', SEPS{SEP_index}, methods{method_index}, N_list{n_data_index}, mean_cost(1), std_cost(1), mean_cost(2), std_cost(2), mean_cost(3), std_cost(3), n_itr_mean, n_itr_std, time_mean, time_std);
            end
        end
    end
    fclose(fid);
end




