%% |positivefactory|
% Returns a manifold struct to optimize over manifold of positive numbers
%
% *Syntax*
%
%   M = positivefactory(n)
%
% *Description*
%
% |M = positivefactory(n)| returns |M|, a structure describing the positive
% manifold.
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Original author: Reshad Hosseini, Dec. 04, 2014.
%
% Change log: 
%

function M = sphere2factory(n)
    
    if ~exist('n', 'var') || isempty(n)
        n = 1;
    end

    M.name = @() sprintf('Simplex space Embedded in R^%i', n);
    
    M.dim = @() n;
    
    M.inner = @(x, d1, d2) d1(:).'*d2(:);
    
    M.norm = @(x, d) norm(d, 'fro');
    
    M.dist = @(x, y) norm(x-y, 'fro');
    
    M.typicaldist = @() sqrt(n);
    
    M.proj = @(x, d) d;
    
    M.egrad2rgrad = @egrad2rgrad;
    function gn = egrad2rgrad(x, g)
        c = x.' * g;
        gn = g - c * x;
    end
    
    %M.ehess2rhess = @(x, eg, eh, d) eh;
    
    M.tangent = M.proj;
    
    M.exp = @expmap;
    function y = expmap(x, d, t)
        % apply first the following change of variable
        %xn = log(x);
        xn = x;
        % moving in the variable change domain
        if nargin == 3
            yn = xn + t*d;
        else
            yn = xn + d;
        end
        % going back to the original domain
        y = yn / norm(yn, 'fro');
    end
    
    M.retr = M.exp;

    M.hash = @(x) ['z' hashmd5(x(:))];
    
    M.rand = @random;
    function x = random()
        x = rand(n,1);
        x = x / norm(x, 'fro');
    end
    
    M.randvec = @randvec;
    function u = randvec(x) %#ok<INUSD>
        u = randn(n,1);
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
