%% MATLAB initials:
clc
clear all
close all

%% installing for adding paths:
addpath(genpath(fullfile("./", 'utils')))
install;

%% settings:
solver_type = "LBFG_VTFree";  %%--> LBFG_manopt_original, LBFG_mixest_original, LBFG_VTFree
manifold_version = "SPD_VTFree";   %%---> SPD_manopt_original, SPD_mixest_original, SPD_mixest_original_fast, SPD_VTFree
experiment = "Karcher_mean";  %%--> Karcher_mean
dimenion_of_matrix = 100;   %%--> 100, 1000, 10000
start_with_given_initial_point = true;
use_saved_initial_point = false;
number_of_runs = 2;
base_dir = "./saved_files/" + experiment + "/n=" + dimenion_of_matrix + "/";

%% optimization runs:
for run_index = 1:number_of_runs
    path_of_initial_point = base_dir + "run" + (run_index) + "/"; 
    path_save = path_of_initial_point + solver_type + "/" + manifold_version + "/";
    if ~exist(path_save, 'dir')
        mkdir(path_save);
    end
    
    %%%%%%%% get optimization problem:
    problem = positive_definite_karcher_mean([], dimenion_of_matrix, manifold_version);
    
    %%%%%%%% generate/get initial point:
    if isfile(path_of_initial_point+"x_initial.mat")  %--> File exists.
        load(path_of_initial_point+"x_initial");
    else   %--> File does not exist.
        x_initial = problem.M.rand();   %--> x_initial can be set to A(:, :, 1) for Karcher_mean experiment
        save(path_of_initial_point+"x_initial.mat", 'x_initial');
    end
    
    %%%%%%%% solve optimization:
    if solver_type == "LBFG_manopt_original"
        if start_with_given_initial_point
            [X, cost_, info_, costevals] = lbfgs_MANOPT(problem, x_initial);
        else
            [X, cost_, info_, costevals] = lbfgs_MANOPT(problem);
        end
    elseif solver_type == "LBFG_mixest_original"
        if start_with_given_initial_point
            [X, cost_, info_, costevals] = lbfgs_MIXEST(problem, x_initial);
        else
            [X, cost_, info_, costevals] = lbfgs_MIXEST(problem);
        end
    elseif solver_type == "LBFG_VTFree"
        if start_with_given_initial_point
            [X, cost_, info_, costevals] = lbfgs_TransportFree(problem, x_initial);
        else
            [X, cost_, info_, costevals] = lbfgs_TransportFree(problem);
        end
    end
    
    %%%%%%%% get the history of optimization:
    [cost_list, grad_norm_list, stepsize_list, time_list, time_iterations] = get_optimization_history(info_);

    %%%%%%%% plot the history of optimization:
    plot_and_save_figure(cost_list, "cost", path_save+"/")
    plot_and_save_figure(log(cost_list), "log of cost", path_save)
    plot_and_save_figure(grad_norm_list, "gradient norm", path_save)
    plot_and_save_figure(log(grad_norm_list), "log of gradient norm", path_save)
    plot_and_save_figure(time_iterations, "time of each itr", path_save)
    plot_and_save_figure(time_list, "time", path_save)
    fprintf('Number of executions of getCostGrad function is : %d\n', costevals);
    save(path_save+'costevals.txt', 'costevals', '-ASCII');

    %%%%%%%% saving the workspace:
    save(path_save+'workspace.mat');
end