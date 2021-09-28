%% |spdfactory - Vector Transport Free|
% Returns a manifold structure to optimize over symmetric positive definite
% matrices
%
% *Syntax*
%
%   M = spdfactory(n)
%
% *Description*
%
% |M = spdfactory(n)| returns |M|, a structure describing the Riemmanian
% manifold of symmetric |n-by-n| positive definite matrices.
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Original author: Suvrit Sra, Aug, 02, 2013
% Contributors: 
%  Reshad Hosseini 
%
% Change log:
%  Reshad Hosseini, Aug,30,2013: Implementing retr, transp, ehess2rhess
%  Reshad Hosseini, Jan,14,2014: Improving speed of transp using sqrtm_fast
%  Reshad Hosseini, Jun,26,2014: Improving mbfgs speed by adding transpF
%  Reza Godaz, Benyamin Ghojogh, April 11, 2021: Changing operations of SPD manifold for vector transport free operations 
%

function M = spdfactory_VTFtreeCholesky(n)   

VTFreeCholesky_flag = true;
global retraction_type;
% retraction_type = "expm"; %--> expm , taylor  ----> we set it in main.m

%%%%%% the flags "flag" and "riemTransp" are ignored if "VTFree_flag" is true
flag = true; % flag = true v. t. riemman ; flag=false: v. t. is identitty
riemTransp = true; % faster to use identity instead of Riemmanian Transp
% If flag is one then it corresponds to transp of natural metrix

if ~exist('n', 'var') || isempty(n)
    n = 1;
end

M.name = @() sprintf('SPD manifold (%d, %d)', n, n);

M.dim = @() (n*(n-1))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--------------------------- basic functions:

M.inner = @inner_product;
    function inner_prod = inner_product(X, U, V)
        if VTFreeCholesky_flag
            inner_prod = U(:).'*V(:);  %--> equaivalent to trace(U.' * V)
        else
            inner_prod = real(sum(sum( (X\U).' .* (X\V) ))); %U(:).'*V(:);
        end
    end

M.norm = @manifold_norm;
    function the_norm = manifold_norm(X, U)
        if VTFreeCholesky_flag
            the_norm = sqrt(U(:).'*U(:));  %--> equaivalent to trace(U.' * V)
        else
            the_norm = sqrt(real(sum(sum( abs(X\U).^2 ))));
        end
    end

M.dist = @riem;
% Riemmanian distance
    function d = riem(X,Y)
        d = eig(X, Y);
        d = norm(log(d));
    end

M.typicaldist = @() sqrt((n*(n-1))/2);

sym = @(X) (X+X')/2;

M.proj = @projection;
    function Up = projection(X, U)
        % Tangent space of symitian matrices is also a symitian matrix
        Up = sym(U);
    end

M.tangent = M.proj;

% For Riemannian submanifolds with euclidean inner product,
% converting a Euclidean gradient into a
% Riemannian gradient amounts to an orthogonal projection.
% Here the inner product is definted as tr(E X^-1 F X^-1). therefore
% We obtain the following for Riemmanian Gradient
M.egrad2rgrad = @egrad2regrad;
    function Up = egrad2regrad(X, U)
        if VTFreeCholesky_flag
            try
                L = chol(X,'lower');
            catch
                epsilon_ = 1e-2;
                L = chol(X + epsilon_ * eye(length(X)),'lower');
            end
            Up = L' * sym(U) * L;
            
        else
            Up = X * sym(U) * X;
            if 0z
                % this gradient corresponding to euclidean innerproduct is slow
                Up = U;
            end
        end
    end

M.ehess2rhess = @ehess2rhess;
    function Hess = ehess2rhess(X, egrad, ehess, eta)
        if VTFreeCholesky_flag
            %TODO
        else
            Hess = X*sym(ehess)*X + 2*sym(H*sym(egrad)*X);
            Hess = Hess - sym(eta*sym(egrad)*X);
        end
    end

M.exp = @exponential;
    function Y = exponential(X, U, t)
        if nargin == 2
            t = 1;
        end
        Y = retraction(X, U, t);
    end

M.log = @logarithm;
    function U = logarithm(X, Y)
        if VTFreeCholesky_flag
            %TODO ---> it might be different
            U = X*logm(X\Y);
            U = sym(U);
        else
            U = X*logm(X\Y);
            U = sym(U);
        end
    end

% M.log = @logarithm;
%     function U = logarithm(X, Y)
%         if VTFreeCholesky_flag
%             %TODO ---> it might be different
%             L = chol(X,'lower');
%             U = X * logm((L')\Y*L');
%             Y = sym(Y);   
%         else
%             U = X*logm(X\Y);
%             U = sym(U);
%         end
%     end

M.hash = @(X) ['z' hashmd5(X(:))];

M.rand = @random;
    function X = random()
        X = randn(n);
        X = (X*X');
    end

M.randvec = @randomvec;
    function U = randomvec(X)
        U = randn(n);
        U = sym(U);
        U = U / norm(U,'fro');
    end

% Linear combination of tangent vectors
M.lincomb = @lincomb;
    function d = lincomb(x, a1, d1, a2, d2) %#ok<INUSL>
        if nargin == 3
            d = a1*d1;
        elseif nargin == 5
            d = a1*d1 + a2*d2;
        else
            error('Bad use of psd.lincomb.');
        end
    end

M.zerovec = @(x) zeros(n);

M.vec = @(x, u_mat) u_mat(:);

M.mat = @(x, u_vec) reshape(u_vec, [n, n]);

M.vecmatareisometries = @() false;

if ~exist('sqrtm_triu_real','file')
    % check if mex files was compiled successfully
    fast_sqrtm = @(x)sqrtm(x);
    warning('sqrtm_triu_real should be compiled to improve performace');
else
    fast_sqrtm = @(x)sqrtm_fast(x);
end

M.map_the_vector = @map_the_vector_;
    function U = map_the_vector_(X, U)
        if VTFreeCholesky_flag
            L=chol(X,'lower');
            U=(L\U)/(L');
            %L_inv = inv(L);
            %U=L_inv*U*L_inv';
        else
            warning('The mode is not vector tranport free but vector mapping is used.');
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--------------------------- retraction and transports:

M.retr = @retraction;
    function Y = retraction(X, U, t)
        if nargin < 3
            t = 1.0;
        end
        if VTFreeCholesky_flag
            L = chol(X,'lower');
            if retraction_type == "expm"
                E = t*U;
                Y = X * expm((L')\E*L');
                %Y = sym(Y);
            elseif retraction_type == "taylor"
                Y = X + L*( t*U + 0.5*t^2* U'*U )*L';
                %Y = sym(Y);
            end
        else
            if flag
                E = t*U;
                Y = X * expm(X\E);
                Y = sym(Y);
            else
                Y = X + t*U;
            end
        end
    end

M.transp = @transpvec;    %--> Benyamin: this is parallel transport in paper 
    function F = transpvec(X, Y, E)
        if VTFreeCholesky_flag
            F = E;
        else
            if flag
                if riemTransp
                    expconstruct= fast_sqrtm(Y/X);  %--> sqrt(Y*inv(X)) --> note: sqrt(Y*inv(X)) is named E in paper
                    F = expconstruct*E*expconstruct';
                else
                    % Identity parallel transport works for LBFGS
                    % There is also proof for the convergence
                    F = E;
                end
            else
                % identity parallel transport
                F = E;
            end
        end
    end

% applying vector transport and save a variable for applying fast version
M.transpstore = @transpvecf;   %--> Benyamin: this function returns sqrt(Y*inv(X)) and inv(sqrt(Y*inv(X))) --> note: sqrt(Y*inv(X)) is named E in paper
    function [expconstruct,iexpconstruct] = transpvecf(X, Y)
        if VTFreeCholesky_flag
            % warning("The function M.transpstore should not be used for VTFree version because it is not useful at all.");
            expconstruct = nan;   %--> no need to it in VTFree version
            iexpconstruct = nan;   %--> no need to it in VTFree version
        else
            if flag
                if riemTransp
                    expconstruct= fast_sqrtm(Y/X);
                    %F = expconstruct*E*expconstruct';
                    if nargout > 1
                       iexpconstruct = inv(expconstruct); 
                    end
                else
                    % Identity parallel transport works for LBFGS
                    % There is also proof for the convergence
                    %F = E;
                    if nargout > 0
                        expconstruct = eye(size(X,1));
                    end
                    if nargout > 1
                        iexpconstruct = eye(size(X,1));
                    end
                end
            else
                % identity parallel transport
                %F = E;
                if nargout > 0
                    expconstruct = eye(size(X,1));
                end
                if nargout > 1
                    iexpconstruct = eye(size(X,1));
                end
            end
        end
    end

% inverse of vector transport
M.itransp = @itranspvec;
    function F = itranspvec(X, Y, E)
        F = transpvec(Y, X, E);
    end

% faster version of vector transport by storing some information
M.transpf = @transpvecfast; 
    function F = transpvecfast(expconstruct, E)
        if VTFreeCholesky_flag
            % warning("The function M.transpf should not be used for VTFree version because it ignores input expconstruct.");
            F = E;
        else
            if flag
                if riemTransp
                    F = expconstruct*E*expconstruct';
                else
                    F = E;
                end
            else
                F = E;
            end
        end
    end
    
% faster version of inverse vector transport by storing some information
M.atranspf = @itranspvecfast; 
    function F = itranspvecfast(iexpconstruct, E)
        if VTFreeCholesky_flag
            % warning("The function M.atranspf should not be used for VTFree version because it ignores input iexpconstruct.");
            F = E;
        else
            if flag
                if riemTransp
                    F = iexpconstruct*E*iexpconstruct';
                else
                    F = E;
                end
            else
                F = E;
            end
        end
    end

end


