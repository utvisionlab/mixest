function M = triufactory(m, n)
% Returns a manifold struct over m-by-n lower-triangular matrices.
%
% function M = triufactory(m, n)
%
% Returns M, a structure describing the Euclidean space of lower-triangular
% m-by-n matrices equipped with the standard Frobenius distance and 
% associated trace inner product as a manifold for Manopt.

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Original author: Reshad Hosseini, Aug, 08, 2015

    
    if ~exist('n', 'var') || isempty(n)
        n = 1;
    end

    M.name = @() sprintf('Euclidean space R^(%dx%d)', m, n);
    
    M.dim = @() m*n;
    
    M.inner = @(x, d1, d2) d1(:).'*d2(:);
    
    M.norm = @(x, d) norm(d, 'fro');
    
    M.dist = @(x, y) norm(x-y, 'fro');
    
    M.typicaldist = @() sqrt(m*n);
    
    M.proj = @(x, d) triu(d);
    
    M.egrad2rgrad = @(x, g) triu(g);
    
    M.ehess2rhess = @(x, eg, eh, d) eh;
    
    M.tangent = M.proj;
    
    M.exp = @exp;
    function y = exp(x, d, t)
        if nargin == 3
            y = x + t*d;
        else
            y = x + d;
        end
        y = triu(y);
    end
    
    M.retr = M.exp;
	
	M.log = @(x, y) triu(y-x);

    M.hash = @(x) ['z' hashmd5(x(:))];
    
    M.rand = @() triu(randn(m, n));
    
    M.randvec = @randvec;
    function u = randvec(x) %#ok<INUSD>
        u = triu(randn(m, n));
        u = u / norm(u, 'fro');
    end
    
    M.lincomb = @lincomb;
    function v = lincomb(x, a1, d1, a2, d2) %#ok<INUSL>
        if nargin == 3
            v = a1*d1;
        elseif nargin == 5
            v = a1*d1 + a2*d2;
        else
            error('Bad usage of euclidean.lincomb');
        end
    end
    
    M.zerovec = @(x) zeros(m, n);
    
    M.transp = @(x1, x2, d) d;
    
    M.pairmean = @(x1, x2) triu(.5*(x1+x2));
    
    M.vec = @(x, u_mat) u_mat(:);
    M.mat = @(x, u_vec) reshape(u_vec, [m, n]);
    M.vecmatareisometries = @() true;

end
