function X = positive_definite_karcher_mean(A)
% Computes a Karcher mean of a collection of positive definite matrices.
%%%%%% for definition of Karcher mean, see --->  https://en.wikipedia.org/wiki/Fr%C3%A9chet_mean

%
% function X = positive_definite_karcher_mean(A)
%
% Input:  A 3D matrix A of size nxnxm such that each slice A(:,:,k) is a
%         positive definite matrix of size nxn.
% 
% Output: A positive definite matrix X of size nxn which is a Karcher mean
%         of the m matrices in A, that is, X minimizes the sum of squared
%         Riemannian distances to the matrices in A:
%            f(X) = sum_k=1^m .5*dist^2(X, A(:, :, k))
%         The distance is defined by the natural metric on the set of
%         positive definite matrices: dist(X,Y) = norm(logm(X\Y), 'fro').
% 
% This simple example is not the best way to compute Karcher means. Its
% purpose it to serve as base code to explore other algorithms. In
% particular, in the presence of large noise, this algorithm seems to not
% be able to reach points with a very small gradient norm. This may be
% caused by insufficient accuracy in the gradient computation.

% This file is part of Manopt and is copyrighted. See the license file.
% 
% Main author: Nicolas Boumal, Sept. 3, 2013
% Contributors:
% 
% Change log:
% 
  
clc
clear all
close all

%%%%%%% settings:
manifold_version = "SPD_mixest_original_fast";   %%---> SPD_manopt_original, SPD_mixest_original, SPD_mixest_original_fast, SPD_VTFree
solver_type = "LBFG_mixest_original";  %%--> LBFG_manopt_original, LBFG_mixest_original, LBFG_VTFree
dimenion_of_matrix = 100;   %%--> 100, 1000, 10000
start_with_given_initial_point = true;
use_saved_initial_point = false;
number_of_runs = 2 ;
base_dir = "./saved_files/n="+dimenion_of_matrix+"/";

    % Generate some random data to test the function if none is given.
    if ~exist('A', 'var') || isempty(A)
%         n = 5;
        n = dimenion_of_matrix;
        m = 50;
        A = zeros(n, n, m);
        ref = diag(max(.1, 1+.1*randn(n, 1)));
        for i = 1 : m
            noise = 0.01*randn(n);
            noise = (noise + noise')/2;
            [V, D] = eig(ref + noise);
            A(:, :, i) = V*diag(max(.01, diag(D)))*V';
        end
    end
    
    % Retrieve the size of the problem:
    % There are m matrices of size nxn to average.
    n = size(A, 1);
    m = size(A, 3);
    assert(n == size(A, 2), ...
           ['The slices of A must be square, i.e., the ' ...
	        'first and second dimensions of A must be equal.']);
    
    % Our search space is the set of positive definite matrices of size n.
    % Notice that this is the only place we specify on which manifold we
    % wish to compute Karcher means. Replacing this factory for another
    % geometry will yield code to compute Karcher means on that other
    % manifold, provided that manifold is equipped with a dist function and
    % a logarithmic map log.
    if manifold_version == "SPD_manopt_original"
        M = sympositivedefinitefactory(n);
    elseif manifold_version == "SPD_mixest_original"
        M = spdffactory(n);
    elseif manifold_version == "SPD_mixest_original_fast"
        M = spdfactory(n);
    elseif manifold_version == "SPD_VTFree"
        M = spdfactory_VTFtree(n);
    end
    
    % Define a problem structure, specifying the manifold M, the cost
    % function and its gradient.
    problem.M = M;
    problem.cost = @cost;
    problem.grad = @grad;
    
    % Explicitly pick an approximate Hessian for the trust-region method.
    % (This is only to show an example of how it can be done; the solver
    % below, rlbfgs, does not use the approximate Hessian; trustregions
    % would.)
%     problem.approxhess = approxhessianFD(problem, struct('stepsize', 1e-4));
    
    % The functions below make many redundant computations. This
    % performance hit can be alleviated by using the caching system. We go
    % for a simple implementation here, as a tutorial example.
    
    % Cost function
    function f = cost(X)
        f = 0;
        for k = 1 : m
            f = f + M.dist(X, A(:, :, k))^2;
        end
        f = f/(2*m);
    end

    % Riemannian gradient of the cost function
    function g = grad(X)
        g = M.zerovec(X);
        for k = 1 : m
            % Update g in a linear combination of the form
            % g = g - [something]/m.
            g = M.lincomb(X, 1, g, -1/m, M.log(X, A(:, :, k)));
        end
    end
    
    function [path_save,x_initial] = make_dir_and_init(base_dir_with_run_no,solver,manifold,problem_manifold)
        path_save= base_dir_with_run_no + solver +"/"+ manifold;
        if ~exist(path_save, 'dir')
            mkdir(path_save);
        end
        if isfile(base_dir_with_run_no+"x_initial.mat")
            % File exists.
            load(base_dir_with_run_no+"x_initial");
        else
             % File does not exist.
            x_initial = problem_manifold.rand();
            save(base_dir_with_run_no+"x_initial.mat", 'x_initial');
        end  
        path_save=path_save+"/";
    end
    %%%%%%%% folder of saving results:

    
    
%     path_save = "./saved_files/n="+n+"/";
%     dircontent = dir(path_save);
%     num_dir = sum([dircontent.isdir]) - 2; %-2 to account for the stupid '.' and '..' returned by dir
%     run_index = num_dir+1;
%     path_save=path_save+"run"+(run_index)+"/";
%     mkdir(path_save);
% 
%     %%%% set initial point:
%     if isfile(path_save+"x_initial.mat")
%          % File exists.
%          load(path_save+"x_initial");
%     else
%          % File does not exist.
%          x_initial = problem.M.rand();
%          save(path_save+"x_initial.mat", 'x_initial');
%     end

    % Execute some checks on the derivatives for early debugging.
    % These things can be commented out of course.
    % The slopes should agree on part of the plot at least. In this case,
    % it is sometimes necessary to inspect the plot visually to make the
    % call, but it is indeed correct.
    % checkgradient(problem);
    % pause;
    
    % Execute this if you want to force using a proper parallel vector
    % transport. This is not necessary. If you omit this, the default
    % vector transport is the identity map, which is (of course) cheaper
    % and seems to perform well in practice.
    % M.transp = M.paralleltransp;
    
    % Issue a call to a solver. Default options are selected.
    % Our initial guess is the first data point. Most solvers work well
    % with this problem. Limited-memory BFGS is one good example:

    for run_index=1:number_of_runs
       base_dir_with_run_no=base_dir+"run"+(run_index)+"/";
       [path_save,x_initial] = make_dir_and_init(base_dir_with_run_no,solver_type,manifold_version,problem.M);
    
        if solver_type == "LBFG_manopt_original"
            if start_with_given_initial_point
                %X = rlbfgs(problem, A(:, :, 1));
                X = rlbfgs(problem, x_initial);
            else
                X = rlbfgs(problem);
            end
        elseif solver_type == "LBFG_mixest_original"
            if start_with_given_initial_point
                [X cost_ info_] = lbfgs_MIXEST(problem, x_initial);
            else
                [X cost_ info_] = lbfgs_MIXEST(problem);
            end
        elseif solver_type == "LBFG_VTFree"
            if start_with_given_initial_point
                [X cost_ info_] = lbfgs_TransportFree(problem, x_initial);
            else
                [X cost_ info_] = lbfgs_TransportFree(problem);
            end
        end


        %%%%%%%% get the history of optimization:
        [cost_list, grad_norm_list, stepsize_list, time_list, time_iterations] = get_optimization_history(info_);

        %%%%%%%% folder of saving results:
        %path_save = "./saved_files/n="+n+"/run"+run_index+"/solver="+solver_type+"/manifold"+manifold_version+"/";
        %if ~exist(path_save, 'dir')
        %    mkdir(path_save);
        %end

        %%%%%%%% plot the history of optimization:
        
        plot_and_save_figure(cost_list, "cost", path_save+"/")
        plot_and_save_figure(log(cost_list), "log of cost", path_save)
        plot_and_save_figure(grad_norm_list, "gradient norm", path_save)
        plot_and_save_figure(log(grad_norm_list), "log of gradient norm", path_save)
        plot_and_save_figure(time_iterations, "time of each itr", path_save)
        plot_and_save_figure(time_list, "time", path_save)

        %%%%%%%% saving the workspace:
        save(path_save+'workspace.mat');
    end
end
