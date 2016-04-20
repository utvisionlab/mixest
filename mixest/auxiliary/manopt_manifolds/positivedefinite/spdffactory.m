%% |spdfastfactory|
% Returns a manifold structure to optimize over symmetric positive definite
% matrices
%
% *Syntax*
%
%   M = spdfastfactory(n)
%
% *Description*
%
% |M = spdfastfactory(n)| returns |M|, a structure describing the 
% Riemmanian manifold of symmetric |n-by-n| positive definite matrices.
%
% The retraction used here is taylor approximation to exponential map
% a function is applied on matrix to make it positive definite
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Original author: Reshad Hosseini, Aug, 02, 2013
%
% Change log:
%

function M = spdffactory(n)   


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
    end

M.ehess2rhess = @ehess2rhess;
    function Hess = ehess2rhess(X, egrad, ehess, eta)
        Hess = X*sym(ehess)*X + 2*sym(H*sym(egrad)*X);
        Hess = Hess - sym(eta*sym(egrad)*X);
    end

% It is possible to apply retraction several times with different t
M.retr = @retraction;
    function [Y, store] = retraction(X, U, t, store)
        if nargin < 3 || isempty(t)
            t = 1.0;
        end     
        if nargin < 4
            store = struct;
        end
        if ~isfield(store, 'UinvXmulU')
            if ~isfield(store, 'invXmulU')
                store.invXmulU = X\U;
            end
            store.UinvXmulU = U*store.invXmulU;
        end
        Y = X + t*U + t^2/2 * store.UinvXmulU;
        Y = sym(Y);        
    end

M.retrtransp = @retractiontranspvec;
    function [Y, F, store] = retractiontranspvec(X, U, E, t, store)
        % retraction at X in direction of U with step-length t
        % vector transport of E into the new point
        if nargin < 4 || isempty(t)
            t = 1.0;
        end   
        if nargin < 5
            store = struct;
        end
        if ~isfield(store, 'UinvXmulU') || ~isfield(store, 'symEinvXmulU')
            if ~isfield(store, 'invXmulU')
                invXmulU = X\U;
            else
                invXmulU = store.invXmulU;
            end
        end
        if ~isfield(store, 'UinvXmulU')
            store.UinvXmulU = U*invXmulU;
        end
        Y = X + t*U + t^2/2 * store.UinvXmulU;
        Y = sym(Y);   
        if ~isfield(store, 'symEinvXmulU')
            store.symUinvXmulU = sym(E * invXmulU);
        end
        F = E + t * store.symUinvXmulU;
    end

% vector transport in line search all arguments except t remains
M.transp = @transpvec;   
    function [F, store] = transpvec(X, Y, E, U, t, store)
        if nargin < 5
            store = struct;
        end
        if ~isfield(store, 'symEinvXmulU')
            if nargin > 3 %~isempty(U)
                if ~isfield(store, 'invXmulU')
                    store.invXmulU = X\U;
                end
                store.symUinvXmulU = sym(E * store.invXmulU);
                F = E + t * store.symUinvXmulU;
            else
                if nargout > 1
                    store.symUinvXmulU = sym((sqrtm(2*Y/X-eye(n))-eye(n))*E);
                    F = E + t * store.symUinvXmulU;
                else
                    F = sym(sqrtm(2*Y/X-eye(n))*E);
                end
            end
        else
            F = E + t * store.symUinvXmulU;
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
    fast_sqrtm = @(x)sqrtm_fast(x);
end

% Applying vector transport aand save a variable for applying fast version
%   of vector transport and its adjoint
% In the fast version only stored variable is given
M.transpstore = @transpvecstore;
    function [F, expconstruct, iexpconstruct] = transpvecstore(X, Y, E)
        expconstruct= fast_sqrtm(Y/X);
        F = expconstruct * E * expconstruct';
        iexpconstruct = inv(expconstruct);
    end

% faster version of vector transport by storing some information
M.transpf = @transpvecfast; 
    function F = transpvecfast(expconstruct, E)
        F = expconstruct * E * expconstruct.';
    end
    
% faster version of adjoint vector transport by storing some information
M.atranspf = @atranspvecfast; 
    function F = atranspvecfast(iexpconstruct, E)
        F = iexpconstruct * E * iexpconstruct.';
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
