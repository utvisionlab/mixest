%% |simplexfactory|
% Returns a manifold struct to optimize over simplex with given dimension
% and sum
%
% *Syntax*
%
%   M = simplexKfactory(n, K)
%
% *Description*
%
% |M = simplexKfactory(n, K)| returns |M|, a structure describing the
% simplex manifold of dimension |n|, where the sum of the elements of the
% points on the manifold equals |K|.
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Original author: Reshad Hosseini, Oct. 29, 2014.
%
% Change log: 
%

function M = simplexKfactory(n, K)
    
    if ~exist('n', 'var') || isempty(n)
        n = 1;
    end

    M.name = @() sprintf('Simplex space Embedded in R^%i', n);
    
    M.dim = @() n-1;
    
    M.inner = @(x, d1, d2) d1(:).'*d2(:);
    
    M.norm = @(x, d) norm(d, 'fro');
    
    M.dist = @(x, y) norm(x-y, 'fro');
    
    M.typicaldist = @() sqrt(n-1);
    
    M.proj = @(x, d) d;
    
    M.egrad2rgrad = @egrad2rgrad;
    function gn = egrad2rgrad(x, g)
        Cov = diag(x) - (x/K) * x';
        gn = Cov(1:end-1,:)*g;
    end
    
    %M.ehess2rhess = @(x, eg, eh, d) eh;
    
    M.tangent = M.proj;
    
    M.exp = @expmap;
    function y = expmap(x, d, t)
        % apply first the following change of variable
        xn = log(x(1:n-1)) - log(x(n));
        % moving in the variable change domain
        if nargin == 3
            yn = xn + t*d;
        else
            yn = xn + d;
        end
        % going back to the original domain
        y = exp([yn;0]);
        y = y / sum(y) * K;
    end
    
    M.retr = M.exp;

    M.hash = @(x) ['z' hashmd5(x(:))];
    
    M.rand = @random;
    function x = random()
        x = rand(n,1);
        x = x / sum(x) * K;
    end
    
    M.randvec = @randvec;
    function u = randvec(x) %#ok<INUSD>
        u = randn(n-1,1);
        u = u / norm(u, 'fro');
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
    
    M.zerovec = @(x) zeros(n, 1);
    
    M.transp = @(x1, x2, d) d;
    
    M.pairmean = @(x1, x2) .5*(x1+x2);
    
    M.vec = @(x, u_mat) u_mat(:);
    M.mat = @(x, u_vec) reshape(u_vec, [m, n]);
    M.vecmatareisometries = @() true;

end
