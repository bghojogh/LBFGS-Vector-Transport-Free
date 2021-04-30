function [X,T,L1,L2,theta] = mixgen(n,m,k,d,c,e)

%mixgen - Gaussian mixture generator
%
%[X,T] = mixgen(n,m,k,d,c,e)
%  n - size of training set
%  m - size of test set
%  k - number of components
%  d - dimension
%  c - separation degree
%  e - maximum eccentricity
%returns
%  X - training set (n x d)
%  T - test set (m x d)

% Nikos Vlassis, 2000
% for definitions see (Dasgupta, 1999)
%
% changed (slightly) by Jan Nunnink, 2003
% changed (slightly) by Reshad Hosseini, 2015

R=zeros(k,d^2);

% mixing weights
while 1
    W = rand(k,1);
    W = W / sum(W);
    if all(W > 1/(4*k))
        break;
    end
end

% create c-separated Gaussian clusters of maximum eccentricity e
trials = 0;
rng('shuffle');
while 1
    X = [];
    T = [];
    M = rand(k,d); % more meaningful to use uniformly sampled data
    %M = bsxfun(@rdivide,M,sqrt(sum(M.^2,2)));
    M = M*sqrt(k)*sqrt(c)*trials/10;
    
    Trace = zeros(k,1);
    Ls = cell(k,1);
    for j = 1:k
        L = rand(d,1)*(e-1)+1;
        if false
            if e ~= 1
                Ls{j} = ( L-min(L) ) / ( max(L)-min(L) ) * (e-1) + 1;
            else
                Ls{j} = L;
            end
        else
            % More meaningful if data exhibits a maximum eccentricity
            % and not all have equal eccentricity
            Ls{j} = L; 
        end
        Trace(j) = sum(Ls{j}.^2/100); 
    end
    
    % check degree of separation
    error = 0;
    for i = 1:k-1
        for j = i+1:k
            if norm(M(i,:)-M(j,:)) < c * sqrt(max(Trace(i),Trace(j)))
                error = 1;
            end
        end
    end
    
    if ~error
        theta.p = W;
        for j = 1:k
            L = diag(Ls{j}).^2/100;
            msg = 1;
            while msg
                % U should be used here so to have different principal axis for covariances
                U = rand(d,d)-0.5;
                U = sqrtm(inv(U*U')) * U;
                [C,msg] = chol(U*L*U');
            end
            Cs = U*L*U';
            R(j,:)=Cs(:)';
            Trace(j) = trace(Cs);
            nj = ceil(n*W(j));
            Xj = randn(nj,d) * C;
            if isempty(X)
                X = repmat(M(j,:),nj,1) + Xj;
            else
                X = [X; repmat(M(j,:),nj,1) + Xj];
            end
            Trace(j) = trace(cov(Xj));
            
            mj = ceil(m*W(j));
            Tj = randn(mj,d) * C;
            if isempty(T)
                T = repmat(M(j,:),mj,1) + Tj;
            else
                T = [T; repmat(M(j,:),mj,1) + Tj];
            end
            %
            % Create theta variable containing initial data
            theta.D{j}.mu = M(j,:);
            theta.D{j}.sigma = Cs;
        end
        break;
    end
    trials = trials + 0.0001;
end


% trials
L = em_gauss(X,M,R);
F = L*W;
F(find(F < eps)) = eps;
L1 = mean(log(F));
if ~isempty(T)
    L = em_gauss(T,M,R);
    F = L*W;
    F(find(F < eps)) = eps;
    L2 = mean(log(F));
end