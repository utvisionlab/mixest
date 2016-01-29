%% |spdfactory|
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
% The retraction used here is taylor approximation to exponential map
% a function is applied on matrix to make it positive definite
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
%

function M = spdfastfactory(n)   

flag = true; % flag = true v. t. riemman ; flag=false: v. t. is identitty
riemTransp = false; % faster to use identity instead of Riemmanian Transp
% If flag is one then it corresponds to transp of natural metrix

if ~exist('n', 'var') || isempty(n)
    n = 1;
end

M.name = @() sprintf('SPD manifold (%d, %d)', n, n);

M.dim = @() (n*(n-1))/2;

M.inner = @(X, U, V) real(sum(sum( (X\U).' .* (X\V) ))); %U(:).'*V(:);

M.norm = @(X, U)  sqrt(real(sum(sum( abs(X\U).^2 ))));

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
        Up = X * sym(U) * X;
        if 0z
            % this gradient corresponding to euclidean innerproduct is slow
            Up = U;
        end
    end

M.ehess2rhess = @ehess2rhess;
    function Hess = ehess2rhess(X, egrad, ehess, eta)
        Hess = X*sym(ehess)*X + 2*sym(H*sym(egrad)*X);
        Hess = Hess - sym(eta*sym(egrad)*X);
    end

M.retr = @retraction;
    function Y = retraction(X, U, t)
        if nargin < 3
            t = 1.0;
        end      
        Y = X + t*U + t^2/2 * U*(X\U);
        Y = sym(Y);        
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
        U = X*logm(X\Y);
        U = sym(U);
    end

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

M.lincomb = @lincomb;

M.zerovec = @(x) zeros(n);

if ~exist('sqrtm_triu_real','file')
    % check if mex files was compiled successfully
    fast_sqrtm = @(x)sqrtm(x);
    warning('sqrtm_triu_real should be compiled to improve performace');
else
    %fast_sqrtm = @(x)sqrtm_fast(x);
    fast_sqrtm = @(x)sqrtm(x);
end
M.transp = @transpvec;   
    function F = transpvec(X, Y, E)
        if flag
            if riemTransp
                expconstruct= fast_sqrtm(Y/X);
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

% applying vector transpord and save a variable for applying fast version
M.transpf = @transpvecf;
    function [F,expconstruct,iexpconstruct] = transpvecf(X, Y, E)
        if flag
            if riemTransp
                expconstruct= fast_sqrtm(Y/X);
                F = expconstruct*E*expconstruct';
                if nargout > 2
                   iexpconstruct = inv(expconstruct); 
                end
            else
                % Identity parallel transport works for LBFGS
                % There is also proof for the convergence
                F = E;
                if nargout > 1
                    expconstruct = eye(size(X,1));
                end
                if nargout > 2
                    iexpconstruct = eye(size(X,1));
                end
            end
        else
            % identity parallel transport
            F = E;
            if nargout > 1
                expconstruct = eye(size(X,1));
            end
            if nargout > 2
                iexpconstruct = eye(size(X,1));
            end
        end
    end

% inverse of vector transport
M.itransp = @itranspvec;
    function F = itranspvec(X, Y, E)
        F = transpvec(Y, X, E);
    end

% faster version of vector transport by storing some information
M.transpF = @transpvecfast; 
    function F = transpvecfast(expconstruct, E)
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
    
% faster version of inverse vector transport by storing some information
M.itranspF = @itranspvecfast; 
    function F = itranspvecfast(iexpconstruct, E)
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

M.vec = @(x, u_mat) u_mat(:);

M.mat = @(x, u_vec) reshape(u_vec, [n, n]);

M.vecmatareisometries = @() false;

end

% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2) %#ok<INUSL>

if nargin == 3
    d = a1*d1;
elseif nargin == 5
    d = a1*d1 + a2*d2;
else
    error('Bad use of psd.lincomb.');
end

end
