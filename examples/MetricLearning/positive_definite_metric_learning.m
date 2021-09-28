function [W, cost_, info_, costevals, ds] = positive_definite_metric_learning(dataset_name, tripletsize_per_class, solver_type, triplet_type, regularization_parameter, pca_n_components)
    %% loaqd dataset: 
    [ds]=load_data(dataset_name);

    %%%%% normalize data:
    if max(max(ds.xTr)) == 255
        ds.xTr=ds.xTr/255.0;
        ds.xTe=ds.xTe/255.0;
    end
    %normalize(ds.xTr,'norm',2);
    %normalize(ds.xTe,'norm',2);
    
    %% apply PCA:
    coeff = pca((ds.xTr)');
    ds.xTr = ds.xTr * coeff(:, 1:pca_n_components);
    ds.xTe = ds.xTe * coeff(:, 1:pca_n_components);
  
    %% create triplets:
    % Generate some random data to test the function if none is given.
    folder_name = pwd + "/examples/MetricLearning/Datasets/triplets/";
    if ~exist(folder_name, 'dir')
       mkdir(folder_name)
    end
    filename_triplets = folder_name + dataset_name + "_" + "TripletsPerClass=" + tripletsize_per_class + '_Triplets.mat';
    if exist(filename_triplets, 'file') == 2
        load(filename_triplets, 'Triplets');
    else
        [Triplets] = CreateTriplets(ds.xTr, ds.yTr, tripletsize_per_class);
        save(filename_triplets, 'Triplets');
    end
    n_triplets = size(Triplets,1);
    X_train = ds.xTr;
    n = size(X_train, 2);

    %% define manifold:
    if solver_type == "VTF_RLBFGS_ISR"
        M = spdfactory_VTFtree(n);
    elseif solver_type == "VTF_RLBFGS_Cholesky"
        M = spdfactory_VTFtreeCholesky(n);
    elseif solver_type == "RLBFGS"
        M = spdfactory_withOptionExpTaylor(n);
    end
    
    %% define cost:
    problem.M = M;
    problem.cost = @cost;
    problem.egrad = @egrad;
    
    % Cost function
    function f = cost(W)
        f = 0;
        for triplet_ind = 1 : n_triplets
            anchor_positive_difference = X_train(Triplets(triplet_ind,1), :) - X_train(Triplets(triplet_ind,2), :);
            f = f + anchor_positive_difference*W*anchor_positive_difference';
            if triplet_type
                anchor_negative_difference = X_train(Triplets(triplet_ind,1), :) - X_train(Triplets(triplet_ind,3), :);
                f = f + anchor_negative_difference * (W \ anchor_negative_difference');
            end
        end
        f = f + 0.5*regularization_parameter*norm(W,'fro');
    end
    
    % Euclidean gradient of the cost function
    function g = egrad(W)
        g = 0;
        for triplet_ind = 1 : n_triplets
            anchor_positive_difference = X_train(Triplets(triplet_ind,1), :) - X_train(Triplets(triplet_ind,2), :);
            g = g + anchor_positive_difference'*anchor_positive_difference;
            if triplet_type
                anchor_negative_difference = X_train(Triplets(triplet_ind,1), :) - X_train(Triplets(triplet_ind,3), :);
                temp = W \ anchor_negative_difference';
                g = g - temp * temp';
            end
        end
        g = g + regularization_parameter*W;
    end

    %% optimization:
    if solver_type == "VTF_RLBFGS_ISR"
        [W, cost_, info_, costevals] = lbfgs_TransportFree(problem);
    elseif solver_type == "VTF_RLBFGS_Cholesky"
        [W, cost_, info_, costevals] = lbfgs_TransportFreeCholesky(problem);
    elseif solver_type == "RLBFGS"
        [W, cost_, info_, costevals] = lbfgs_MIXEST(problem);
    end
    
end
