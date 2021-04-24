function L = em_gauss(X,M,R)

% em_gauss - compute likelihoods for all points and all components
%
% L = em_gauss(X,M,R)
%  X - (n x d) matrix of input data
%  M - (k x d) matrix of components means
%  R - (k x d^2) matrix of components covariances in vector reshaped format.
% returns 
%  L - (n x k) likelihoods of points x_n belonging to component k
%
% Jan Nunnink, 2003

[n,d] = size(X);
k = size(M,1);

L = zeros(n,k); 
for j = 1:k 

  % Cholesky triangular matrix of component's covariance matrix
  Rj = reshape(R(j,:),d,d);
  Rj = chol(Rj);
  
  % We need to compute the Mahalanobis distances between all inputs
  % and the mean of component j; using the Cholesky form of covariances
  % this becomes the Euclidean norm of some new vectors 
  New = (X - repmat(M(j,:),n,1)) * inv(Rj);
  Mah = sum2(New.^2,2);

  L(:,j) = (2*pi)^(-d/2) / det(Rj) * exp(-0.5*Mah);
end

function S = sum2(A,dim);

% calculate sum of A over dimension dim
%
% Jan Nunnink, 2002

[d1, d2] = size(A);

if dim==1
   if d1==1
      S=A;
   else
      S=sum(A);
   end
else
   if d2==1
      S=A;
   else
      S=sum(A')';
   end
end

