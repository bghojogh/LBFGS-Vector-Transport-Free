function theta0 = init_param(data, K, INIT)
D = mvnfactory(size(data,1));
D2 = mixturefactory(D, K);
options = mxe_options();
options.penalize = true;
options.penalizertheta = D2.penalizerparam(data);
LL_save = -inf;
if K > 1
    iter = 1;
    while true
        switch INIT
            case 'default'
                theta0 = [];
                return
                
            case 'kmeans'
                [L, C] = kmeans(data.', K);
                L = L.'; C = C.';
                
            case 'kmeanspp'
                [L, C] = kmeanspp(data, K);
        end
        
        % convert k-means output to mixture params
        theta0.p = zeros(K,1);
        for k = 1:K
            theta0.D{k}.mu = C(:,k);
            X = data(:,L==k);
            X = X - repmat(C(:,k), [1,size(X,2)]);
            theta0.D{k}.sigma = X * X.' / size(X,2);
            theta0.p(k) = sum(L==k) / length(L);
            [R, p] = chol(theta0.D{k}.sigma);
            lenComponent = sum(L==k);
            lenDataPerC = size(data,2) / K;
            if lenComponent < 0.1 * lenDataPerC
                disp('Component have small percentage of data');
                p = 1;
                break;
            end
            if lenComponent < size(data,1)
                disp('Length is less that number of data points');
                p = 1;
                break;
            end
            if p ~= 0
                disp('not spd');
                break;
            end
        end
        
        if p == 0
            iter = iter + 1;
            cost = mxe_costgrad( D2, theta0, data, options );
            LL = -1 * cost;
            if LL > LL_save && LL ~= inf
                LL_save = LL;
                theta0Save = theta0;
            end
        end
        if iter > 30
            break;
        end
    end
    theta0 = theta0Save;
else
    theta0.mu = data(:,1);
    data = data - repmat(theta0.mu, [1,size(data,2)]);
    theta0.sigma = data * data.' / size(data,2);
end
