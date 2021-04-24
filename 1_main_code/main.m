%% MATLAB initials:
clc
clear all
close all

%% installing for adding paths:
addpath(genpath(fullfile("./", 'utils')))
install;

%% settings:
solver_type = "RLBFGS_Wolfe";  %%--> RLBFGS_cautious, RLBFGS_Wolfe, RLBFGS_Wolfe_VTFree, RLBFGS_Wolfe_VTFreeCholesky
global retraction_type; retraction_type = "taylor"; %--> expm , taylor
experiment = "Karcher_mean";  %%--> Karcher_mean
dimenion_of_matrix = 100;   %%--> 100, 1000, 10000
start_with_given_initial_point = true;
number_of_runs = 2;

%% set the manifold based on the solver:
%%---> SPD_manopt_original, SPD_mixest_original, SPD_mixest_original_fast, SPD_VTFree, SPD_VTFreeCholesky
if solver_type == "RLBFGS_cautious"
    manifold_version = "SPD_manopt_original"; 
elseif solver_type == "RLBFGS_Wolfe" && retraction_type == "taylor"
    manifold_version = "SPD_mixest_original"; 
elseif solver_type == "RLBFGS_Wolfe" && retraction_type == "expm"
    manifold_version = "SPD_mixest_original_fast"; 
elseif solver_type == "RLBFGS_Wolfe_VTFree"
    manifold_version = "SPD_VTFree"; 
elseif solver_type == "RLBFGS_Wolfe_VTFreeCholesky"
    manifold_version = "SPD_VTFreeCholesky"; 
end

%% optimization runs:
base_dir = "./saved_files/" + experiment + "/n=" + dimenion_of_matrix + "/";
for run_index = 1:number_of_runs
    path_of_initial_point = base_dir + "run" + (run_index) + "/"; 
    if solver_type == "RLBFGS_Wolfe_VTFree" || solver_type == "RLBFGS_Wolfe_VTFreeCholesky" || solver_type == "RLBFGS_Wolfe"
        solver_type_ = solver_type + "_" + retraction_type;
    else
        solver_type_ = solver_type;
    end
    path_save = path_of_initial_point + solver_type_ + "/" + manifold_version + "/";
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
    if solver_type == "RLBFGS_cautious"   %--> LBFGS_manopt_original
        if start_with_given_initial_point
            [X, cost_, info_, op , costevals] = lbfgs_MANOPT(problem, x_initial);
        else
            [X, cost_, info_, op , costevals] = lbfgs_MANOPT(problem);
        end
    elseif solver_type == "RLBFGS_Wolfe"   %--> LBFGS_mixest_original
        if start_with_given_initial_point
            [X, cost_, info_, costevals] = lbfgs_MIXEST(problem, x_initial);
        else
            [X, cost_, info_, costevals] = lbfgs_MIXEST(problem);
        end
    elseif solver_type == "RLBFGS_Wolfe_VTFree"
        if start_with_given_initial_point
            [X, cost_, info_, costevals] = lbfgs_TransportFree(problem, x_initial);
        else
            [X, cost_, info_, costevals] = lbfgs_TransportFree(problem);
        end
    elseif solver_type == "RLBFGS_Wolfe_VTFreeCholesky"
        if start_with_given_initial_point
            [X, cost_, info_, costevals] = lbfgs_TransportFreeCholesky(problem, x_initial);
        else
            [X, cost_, info_, costevals] = lbfgs_TransportFreeCholesky(problem);
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