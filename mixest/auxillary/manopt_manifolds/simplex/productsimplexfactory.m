%% |simplexfactory|
% Returns a manifold struct to optimize over simplex with given dimension
%
% *Syntax*
%
%   M = simplexfactory(n)
%
% *Description*
%
% |M = simplexfactory(n)| returns |M|, a structure describing the simplex
% manifold of dimension |n|.
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Original author: Reshad Hosseini, Oct. 29, 2014.
%
% Change log: 
%

function M = productsimplexfactory(n, m)

    if nargin < 2
        m = 1;
    end
    
    if ~exist('n', 'var') || isempty(n)
        n = 1;
    end

    M.name = @() sprintf('Simplex space Embedded in R^%i', n);
    
    M.dim = @() (n-1)* m ;
    
    M.inner = @(x, d1, d2) d1(:).'*d2(:);
    
    M.norm = @(x, d) norm(d, 'fro');
    
    M.dist = @(x, y) norm(x-y, 'fro');
    
    M.typicaldist = @() sqrt(n-1);
    
    M.proj = @(x, d) d;
    
    M.egrad2rgrad = @egrad2rgrad;
    function gn = egrad2rgrad(x, g)
        gn = zeros(n-1, m);
        for k = 1:m
            Cov = diag(x(:,k)) - x(:,k) * x(:,k)';
            gn(:,k) = Cov(1:end-1,:)*g(:,k);
        end
    end
    
    %M.ehess2rhess = @(x, eg, eh, d) eh;
    
    M.tangent = M.proj;
    
    M.exp = @expmap;
    function y = expmap(x, d, t)
        if size(x, 1) == 1
            y = ones(n, m);
            return;
        end
        % apply first the following change of variable
        xn = bsxfun(@minus, log(x(1:n-1, :)), log(x(n, :)));
        % moving in the variable change domain
        if nargin == 3
            yn = xn + t*d;
        else
            yn = xn + d;
        end
        % going back to the original domain
        y = exp([yn;zeros(1,m)]);
        y = bsxfun(@rdivide, y, sum(y));
    end
    
    M.retr = M.exp;

    M.hash = @(x) ['z' hashmd5(x(:))];
    
    M.rand = @random;
    function x = random()
        x = rand(n,m);
        x = bsxfun(@rdivide, x, sum(x));
    end
    
    M.randvec = @randvec;
    function u = randvec(x) %#ok<INUSD>
        u = randn(n-1,m);
        u = bsxfun( @rdivide, u, sqrt(sum(u.^2)) );
    end
    
    M.lincomb = @lincomb;
    function v = lincomb(x, a1, d1, a2, d2) %#ok<INUSL>
        if nargin == 3
            v = a1*d1;
        elseif nargin == 5
            v = a1*d1 + a2*d2;
        else
            error('Bad usage of simplex.lincomb');
        end
    end
    
    M.zerovec = @(x) zeros(n-1, 1);
    
    M.transp = @(x1, x2, d) d;
    
    M.pairmean = @(x1, x2) .5*(x1+x2);
    
    M.vec = @(x, u_mat) u_mat(:);
    M.mat = @(x, u_vec) reshape(u_vec, [m, n]);
    M.vecmatareisometries = @() true;

end
