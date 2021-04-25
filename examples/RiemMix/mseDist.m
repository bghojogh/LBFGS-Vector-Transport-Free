% function for calcularing MSE distance between mean and covariance matrices
function [mu_dists, sigma_dists] = mseDist(theta1, theta2)
K = length(theta1.D);
if isfield (theta1.D{1},'sigmat')
    D = mvn2factory(size(theta1.D{1}.sigmat,1)-1);
    for k1 = 1:K
        theta1.D{k1} = D.new2main(theta1.D{k1});
    end
end
if isfield (theta2.D{1},'sigmat')
    D = mvn2factory(size(theta2.D{1}.sigmat,1)-1);
    for k1 = 1:K
        theta2.D{k1} = D.new2main(theta2.D{k1});
    end
end

% print mus and sigmas
for k1 = 1:K
    theta1.D{k1}.mu;
    theta1.D{k1}.sigma;
    theta2.D{k1}.mu;
    theta2.D{k1}.sigma;
end
mu_dists = zeros(K);
sigma_dists = zeros(K);
for k1 = 1:K
    for k2 = 1:K
        mu_dists(k1,k2) = norm(theta1.D{k1}.mu - theta2.D{k2}.mu.', 'fro');
        sigma_dists(k1,k2) = norm(theta1.D{k1}.sigma - theta2.D{k2}.sigma, 'fro');
    end
end