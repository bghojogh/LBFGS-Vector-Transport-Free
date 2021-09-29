%% MATLAB initials:
clc
clear all
close all

%% installing for adding paths:
addpath(genpath(fullfile("./", 'utils')))
install;

%% Settings:
experiment = "MetricLearning";  
dataset_name = 'fisheriris';  %--> usps, vehicle, mnist, fisheriris
print_table_again = true;
average_results_again = true;

%% load results
base_dir = "./saved_files/" + experiment + "/";
path = sprintf("%sdataset=%s/", base_dir, dataset_name)+"/all_info.mat";
load(path);  %--> load all_info.mat

%% average results
solver_type_list = ["VTF_RLBFGS_ISR" , "VTF_RLBFGS_Cholesky" , "RLBFGS"];
n_runs = 10;
iterations_to_report = {"last"};  %--> {10, 20, "last"}, {20, 50, "last"}, {30, 50, "last"}, {"last"}--> if using this, we should the headers and computations of times in fprintf and below code
if average_results_again
    costs_list = zeros(n_runs, length(solver_type_list), length(iterations_to_report));
    n_iterations_list = zeros(n_runs, length(solver_type_list));
    time_list = zeros(n_runs, length(solver_type_list));
    time_average_list = zeros(n_runs, length(solver_type_list));
    for run = 1: n_runs
        for solver_ind = 1:length(solver_type_list)
            info_experiment = all_info(run,solver_ind);
            cost = convert_struct_to_array(info_experiment.info_, "cost");
            time = convert_struct_to_array(info_experiment.info_, "time");
            for itr_of_cost_index = 1:length(iterations_to_report)
                itr_ = iterations_to_report{itr_of_cost_index};
                if isstring(itr_)
                    costs_list(run, solver_ind, itr_of_cost_index) = cost(end);
                else
                    if itr_ > length(cost); continue; end
                    costs_list(run, solver_ind, itr_of_cost_index) = cost(itr_);
                end
            end
            n_iterations_list(run, solver_ind) = length(cost);
            time_list(run, solver_ind) = sum(time);
            time_average_list(run, solver_ind) = mean(time);
        end
    end
    
    %--> take mean and std:
    if print_table_again
        path_ = sprintf("%sdataset=%s/", base_dir, dataset_name);
        fid = fopen(path_+'/results.txt', 'wt');
        fprintf(fid, 'solver \t n_iters \t time \t average time(per itr) \t cost(itr=last) \n');
        fprintf(fid, '=============================================== \n');
        %--> save results in table format in text file:
        for solver_ind = 1:length(solver_type_list)
            for itr_of_cost_index = 1:length(iterations_to_report)
                mean_cost(itr_of_cost_index) = mean(costs_list(:, solver_ind, itr_of_cost_index));
                std_cost(itr_of_cost_index) = std(costs_list(:, solver_ind, itr_of_cost_index));
            end
            n_itr_mean = mean(n_iterations_list(:, solver_ind));
            n_itr_std = std(n_iterations_list(:, solver_ind));
            time_mean = mean(time_list(:, solver_ind));
            time_std = std(time_list(:, solver_ind));
            time_average_mean = mean(time_average_list(:, solver_ind));
            time_average_std = std(time_average_list(:, solver_ind));
            %fprintf(fid, '%s \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \n', solver_type_list(solver_ind), mean_cost(1), std_cost(1), mean_cost(2), std_cost(2), mean_cost(3), std_cost(3), n_itr_mean, n_itr_std, time_mean, time_std, time_average_mean, time_average_std);
            fprintf(fid, '%s \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \t %.3f|+|%.3f \n', solver_type_list(solver_ind), n_itr_mean, n_itr_std, time_mean, time_std, time_average_mean, time_average_std, mean_cost(1), std_cost(1));
        end
        fclose(fid);
    end
    
    %%%%%%%%%%%%%%%%%%%%% table for time difference:
    path_ = sprintf("%sdataset=%s/", base_dir, dataset_name);
    fid = fopen(path_+'/results_timeDiff.txt', 'wt');
%     fprintf(fid, 'Method \t %s \t %s \n', solver_type_list(2), solver_type_list(3));  %---> method indices are based on the solver_type_list above
    fprintf(fid, 'Method \t %s \n', solver_type_list(3));  %---> method indices are based on the solver_type_list above
    fprintf(fid, '=============================================== \n');
    %--> save results in table format in text file:
    for method_index1 = 1:length(solver_type_list)
        if ~strcmp(solver_type_list(method_index1), 'VTF_RLBFGS_ISR') && ~strcmp(solver_type_list(method_index1), 'VTF_RLBFGS_Cholesky')
            continue
        end
        fprintf(fid, '\t %s', solver_type_list(method_index1));
        for method_index2 = 1:length(solver_type_list)
            if strcmp(solver_type_list(method_index2), 'VTF_RLBFGS_ISR') || strcmp(solver_type_list(method_index2), 'VTF_RLBFGS_Cholesky')
                continue
            end
            time_diff_mean = mean(time_list(:, method_index1) - time_list(:, method_index2));
            time_diff_std = std(time_list(:, method_index1) - time_list(:, method_index2));
            % fprintf(fid, '\t %.3f|+|%.3f', time_diff_mean, time_diff_std);
            fprintf(fid, '\t %.3f', time_diff_std);
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end




