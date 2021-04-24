function sim1_gen_table(iter_End, DIM, SEPS, K, NDIM, E, INIT, SELECT)

addpath(genpath('/thirdparty/'))
clc
close all

warning off all

cs   = [  0  ,  1  ,  5   ];


% We have only one run for each initialization
run = 1;


% optimization methods
%METHODS = get_methods_verbeek(DIM, K);
%methods = fieldnames(METHODS);

%methods = {'EM', 'SGD', 'LBFGS', 'CG'};
methods = {'LBFGS'};

fprintf('============\n')

for isep = 1:numel(SEPS)
    SEP = SEPS{isep};
    c = cs(isep);
    
    fprintf('\\rowcolor{Gray}\n');
    %fprintf('$c=%d$ ', c)
    fprintf('%s', SEP)
    
    for imethod = 1:numel(methods)
        if imethod > 1
            fprintf('\\rowcolor{Gray}\n');
        end
        
        methodi = methods{imethod};
        fprintf('& %s ', methodi);
        for N = (NDIM* DIM^2)
%             if strcmp(methodi, 'EM')
%                 submethods = {'EM1'};
%             elseif strcmp(methodi, 'LBFGS')
%                 submethods = {'LBFGS2'};
%             elseif strcmp(methodi, 'CG')
%                 submethods = {'CG2'};
%             else
%                 submethods = {'SGDf11', 'SGDf9', 'SGDf12'};
%             end
            %if strcmp(methodi, 'EM')
            %    submethods = {'EM1'};
            if strcmp(methodi, 'LBFGS')
                submethods = {'LBFGS1','LBFGS2'};
            %elseif strcmp(methodi, 'CG')
            %    submethods = {'CG2'};
            %else
            %    submethods = {'SGDf11', 'SGDf9', 'SGDf12'};
            end

            for isubmethod = 1:numel(submethods)
                method = submethods{isubmethod};
                
                if strcmp(method, 'SGDf12') || strcmp(method, 'EM1') || strcmp(method, 'LBFGS2') || strcmp(method, 'CG2')
                    if N > NDIM(1) * DIM^2
                        fprintf('& ');
                    end
                    fprintf('& $%d$ ', N)
                end
                if strcmp(method, 'SGDf11')
                    mses = zeros(3, iter_End);
                end
                
                if strcmp(method, 'EM1') || strcmp(method, 'LBFGS2') || strcmp(method, 'CG2')
                    mses = zeros(4, iter_End);
                    
                end
                
                for iter = 1 : iter_End
                    
                    RESFOLDER = ['result' num2str(iter)];
                    
                    [time_run, iter_run, pll_run, ll_run, mu_mse_run, sigma_mse_run] = ...
                        sim1_cost_mse_likelihood_results(method, DIM, SEP, INIT, K, N, E, run, RESFOLDER);
                    
                    % change sigma to mu
                    if strcmp(SELECT, 'PLL')
                        mu_mse_run = pll_run;
                    elseif strcmp(SELECT, 'SIGMA');
                        mu_mse_run = sigma_mse_run;
                    end
                    
                    if strcmp(method, 'SGDf11')
                        mses(3, iter) = mu_mse_run(end);
                    end
                    if strcmp(method, 'SGDf9')
                        mses(2, iter) = mu_mse_run(end);
                    end
                    if strcmp(method, 'SGDf12')
                        mses(1, iter) = mu_mse_run(end);
                    end
                    if strcmp(method, 'EM1') || strcmp(method, 'LBFGS2') || strcmp(method, 'CG2')
                        
                        mses(1:3, iter) = interp1(iter_run, mu_mse_run, [5 20 50]);
                        if iter_run(end) < 5
                            mses(1, iter) = mu_mse_run(end);
                            mses(2, iter) = mu_mse_run(end); %inf;
                            mses(3, iter) = mu_mse_run(end); %inf;
                        elseif iter_run(end) < 20
                            mses(3, iter) = mu_mse_run(end);% inf;
                            mses(2, iter) = mu_mse_run(end);
                        elseif iter_run(end) < 50
                            mses(3, iter) = mu_mse_run(end);
                        end
                        if iter_run(2) > 5
                            mses(1, iter) = mu_mse_run(2);
                        end
                        mses(4, iter) = mu_mse_run(end);
                        
                    end
                    
                end
            end
            if strcmp(method, 'SGDf12')
                if strcmp(SELECT, 'PLL')
                    fprintf('& %.3f & %.3f & %.3f ', ...
                        -mean(mses(1,:)), -mean(mses(2,:)), -mean(mses(3,:)));
                else
                    fprintf('& %.3f$\\pm$%.3f & %.3f$\\pm$%.3f & %.3f$\\pm$%.3f' , ...
                        mean(mses(1,:)), std(mses(1,:)), mean(mses(2,:)), std(mses(2,:)), mean(mses(3,:)), std(mses(3,:)));
                end
                fprintf('\\\\\n')
            end
            if strcmp(method, 'EM1') || strcmp(method, 'LBFGS2') || strcmp(method, 'CG2')
                if strcmp(SELECT, 'PLL')
                    fprintf('& %.3f & %.3f & %.3f & %.3f', ...
                        -mean(mses(1,:)), -mean(mses(2,:)), -mean(mses(3,:)), -mean(mses(4,:)))
                else
                    fprintf('& %.3f$\\pm$%.3f & %.3f$\\pm$%.3f & %.3f$\\pm$%.3f & %.3f$\\pm$%.3f' , ...
                        mean(mses(1,:)), std(mses(1,:)), mean(mses(2,:)), std(mses(2,:)), mean(mses(3,:)), std(mses(3,:)), mean(mses(4,:)), std(mses(4,:)))
                end
                fprintf('\\\\\n')
            end
            
        end
        
    end
    fprintf('\\hline \n')
end



