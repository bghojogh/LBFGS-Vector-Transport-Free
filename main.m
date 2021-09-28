%% MATLAB initials:
clc
clear all
close all

%% installing for adding paths:
addpath(genpath(fullfile("./", 'utils')))
install;

%% general settings:
global retraction_type; retraction_type = "expm"; %--> expm, taylor
experiment = "MetricLearning";  %%--> Karcher_mean, RiemMix , MetricLearning
number_of_runs = 10;

%% settings for Karcher mean:
if experiment == "Karcher_mean"
    dimenion_of_matrix = 100;
    solver_type = "RLBFGS_Wolfe_VTFree";  %%--> RLBFGS_cautious, RLBFGS_Wolfe, RLBFGS_Wolfe_VTFree, RLBFGS_Wolfe_VTFreeCholesky
    start_with_given_initial_point = true;
    %%%%% set the manifold based on the solver: ---> SPD_manopt_original, SPD_mixest_original, SPD_mixest_original_fast, SPD_VTFree, SPD_VTFreeCholesky
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
end

%% settings for RiemMix: 
%%%%--> note: For the RiemMix experiment, code solves using all algorithms (solvers) (by for loops) and saves the results
if experiment == "RiemMix"
    dimenion_of_matrix = 2;   %%--> for RiemMix: 2 / for Karcher_mean: 100
    DIMS = [dimenion_of_matrix]; % Dimension
    SEPS = {'low','mid','high'}; % Separation
    KS = [2]; % Number of Components
    NDIM = [100]; % Number of Data = NDIM*(DIM^2)  ---> example: [10 100]
    ES = [10]; % Eccentricity
    INITS = {'kmeanspp'}; % Initialization --> 'default', 'kmeans', 'kmeanspp'
    SELECT = 'PLL'; % 'SIGMA', 'PLL', 'MU'
end

%% settings for MetricLearning: 
if experiment == "MetricLearning"
    dataset_name = 'iris';  %--> usps, vehicle, mnist, iris
    tripletsize_per_class = 20;
    solver_type = "VTF_RLBFGS_ISR";  % VTF_RLBFGS_ISR , VTF_RLBFGS_Cholesky , RLBFGS
    triplet_type = true;
    regularization_parameter = 1;
    pca_n_components = 10; 
end

%% optimization runs:
base_dir = "./saved_files/" + experiment + "/";
info_ind = 1;
for run_index = 1:number_of_runs
    if experiment == "Karcher_mean"
        path_of_initial_point = base_dir + "dim=" + dimenion_of_matrix + "/run" + (run_index) + "/"; 
        if solver_type == "RLBFGS_Wolfe_VTFree" || solver_type == "RLBFGS_Wolfe_VTFreeCholesky" || solver_type == "RLBFGS_Wolfe"
            solver_type_ = solver_type + "_" + retraction_type;
        else
            solver_type_ = solver_type;
        end
        path_save = path_of_initial_point + solver_type_ + "/" + manifold_version + "/";
        if ~exist(path_save, 'dir')
            mkdir(path_save);
        end
        %%%%%%%% recording command window:
        record_command_window(path_save, "on")
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
        plot_results("", info_, path_save, costevals)
        record_command_window(path_save, "off")
        all_info(run_index) = {info_};
    elseif experiment == "RiemMix"
        path_save = sprintf("%sdim=%d/run%d/", base_dir, dimenion_of_matrix, run_index);
        if ~exist(path_save, 'dir')
            mkdir(path_save);
        end
        record_command_window(path_save+"plots2/", "on")
        info_list = sim1(1, DIMS, SEPS, KS, NDIM, ES, INITS, SELECT, path_save);
        for i = info_ind:size(info_list,2)+info_ind-1  %length(fieldnames(info_list))-1
            temp_ind = i-info_ind+1;
            run = "run" + run_index;
            info_list(temp_ind).run = run;
            all_info(i) = info_list(temp_ind);
            all_info(i).run= run;
            name = info_list(temp_ind).separation;
            n_data = info_list(temp_ind).n_data;
            select_mode = info_list(temp_ind).select_mode;
            sub_name = "sim1_init_dim("+info_list(temp_ind).dim +")_N("+n_data+")_sep("+name+")_mode("+select_mode+")";
            for j=1: size(info_list(temp_ind).info_list,2)
                plot_results(sub_name,info_list(temp_ind).info_list(j).info, path_save+"/plots2/"+info_list(temp_ind).info_list(j).method+"/");
            end
        end
        info_ind = info_ind + size(info_list,2);
        record_command_window(path_save+"plots2/", "off");

    elseif experiment == "MetricLearning"
        path_save = sprintf("%sdim=%d/run%d/", base_dir, dataset_name);
        if ~exist(path_save, 'dir')
            mkdir(path_save);
        end
        record_command_window(path_save+"plots2/", "on");
        [W, cost_, info_, costevals, ds] = positive_definite_metric_learning(dataset_name, tripletsize_per_class, solver_type, triplet_type, regularization_parameter, pca_n_components);
        %%%% KNN before metric learning:
        Mdl = fitcknn(ds.xTr, ds.yTr, 'NumNeighbors',10);
        [label_pred_tr1,score,cost] = predict(Mdl, ds.xTr);
        accuracy_train1 = sum(label_pred_tr1 == ds.yTr) / length(ds.yTr);
        [label_pred_te1,score,cost] = predict(Mdl, ds.xTe);
        accuracy_test1 = sum(label_pred_te1 == ds.yTe) / length(ds.yTe);
        disp("Train accuracy (before metric learning): " + accuracy_train1);
        disp("Test accuracy (before metric learning): " + accuracy_test1);
        %%%% KNN after metric learning:
        X_train_projected_rowWise = (W' * (ds.xTr)')';
        X_test_projected_rowWise = (W' * (ds.xTe)')';
        Mdl = fitcknn(X_train_projected_rowWise, ds.yTr, 'NumNeighbors',10);
        [label_pred_tr2,score,cost] = predict(Mdl, X_train_projected_rowWise);
        accuracy_train2 = sum(label_pred_tr2 == ds.yTr) / length(ds.yTr);
        [label_pred_te2,score,cost] = predict(Mdl, X_test_projected_rowWise);
        accuracy_test2 = sum(label_pred_te2 == ds.yTe) / length(ds.yTe);
        disp("Train accuracy (after metric learning): " + accuracy_train2);
        disp("Test accuracy (after metric learning): " + accuracy_test2);
        record_command_window(path_save+"plots2/", "off");
    end
    %%%%%%%% saving the workspace:
    save(path_save+"workspace.mat");
end
if experiment == "Karcher_mean"
    save(pwd+"/saved_files/"+experiment+"/dim="+dimenion_of_matrix+"/all_info_"+solver_type_+".mat","all_info");
elseif experiment == "RiemMix"
    save(pwd+"/saved_files/"+experiment+"/all_info.mat","all_info");
end


%% some functions:
function plot_results(name, info_, path_save, costevals)
    %%%%%%%% get the history of optimization:
    [cost_list, grad_norm_list, stepsize_list, time_list, time_iterations] = get_optimization_history(info_);
    %%%%%%%% plot the history of optimization:
    plot_and_save_figure(cost_list, "cost_"+name, "cost", path_save);
    plot_and_save_figure(log(cost_list), "log of cost_"+name, "log of cost", path_save);
    plot_and_save_figure(grad_norm_list, "gradient norm_"+name, "gradient norm", path_save);
    plot_and_save_figure(log(grad_norm_list), "log of gradient norm_"+name, "log of gradient norm", path_save);
    plot_and_save_figure(time_iterations, "time of each itr_"+name, "time of each iteration", path_save);
    plot_and_save_figure(time_list, "time_"+name, "time", path_save);
    plot_and_save_figure(stepsize_list, "step_size_"+name, "step size", path_save);
    if nargin > 3
        fprintf('Number of executions of getCostGrad function is : %d\n', costevals);
        save(path_save+'costevals.txt', 'costevals', '-ASCII');
    end
end

