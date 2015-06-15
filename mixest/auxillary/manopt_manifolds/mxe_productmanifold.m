%% |mxe_productmanifold|
% This function presents an alternative to the |productmanifold| function
% in the Manopt toolbox, using cell arrays instead of structures.
%
% *Syntax*
%
%   M = mxe_productmanifold(elements)
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Original author: Reshad Hosseini
%
% Change log: 
%

function M = mxe_productmanifold(elements)
% This function presents an alternative to the productmanifold function in
% the Manopt toolbox, using cell arrays instead of structures.
%
% Input:
% elements: a cell array of manifold structures.
%
% Output:
% M: a manifold structure representing the Cartesian product of the
% manifolds given in elements, with points and vectors represented as cell
% arrays.


    %elems = fieldnames(elements);
    nelems = numel(elements);

    assert(nelems >= 1, ...
           'elements must be a structure with at least one field.');
    
    M.name = @name;
    function str = name()
        str = 'Product manifold: [';
        
        for i = 1 : nelems-1
            str = [str sprintf(' %s x', ...
                   elements{i}.name())]; %#ok<AGROW>
        end
        str = [str sprintf(' %s ] ', ...
                   elements{end}.name())];
    end
    
    M.dim = @dim;
    function d = dim()
        d = 0;
        for i = 1 : nelems
            d = d + elements{i}.dim();
        end
    end
    
    M.inner = @inner;
    function val = inner(x, u, v)
        val = 0;
        for i = 1 : nelems
            val = val + elements{i}.inner(x{i}, ...
                                               u{i}, v{i});
        end
    end

    M.norm = @(x, d) sqrt(M.inner(x, d, d));

    M.dist = @dist;
    function d = dist(x, y)
        sqd = 0;
        for i = 1 : nelems
            sqd = sqd + elements{i}.dist(x{i}, ...
                                                 y{i})^2;
        end
        d = sqrt(sqd);
    end
    
    M.typicaldist = @typicaldist;
    function d = typicaldist
        sqd = 0;
        for i = 1 : nelems
            sqd = sqd + elements{i}.typicaldist()^2;
        end
        d = sqrt(sqd);
    end

    M.proj = @proj;
    function v = proj(x, u)
        v = cell(nelems, 1);
        for i = 1 : nelems
            v{i} = elements{i}.proj(x{i}, ...
                                                    u{i});
        end
    end

    M.tangent = @tangent;
    function v = tangent(x, u)
        v = cell(nelems, 1);
        for i = 1 : nelems
            v{i} = elements{i}.tangent(x{i}, ...
                                                       u{i});
        end
    end

    M.tangent2ambient = @tangent2ambient;
    function v = tangent2ambient(x, u)
        v = cell(nelems, 1);
        for i = 1 : nelems
            if isfield(elements{i}, 'tangent2ambient')
                v{i} = ...
                    elements{i}.tangent2ambient( ...
                                               x{i}, u{i});
            else
                v{i} = u{i};
            end
        end
    end

    M.egrad2rgrad = @egrad2rgrad;
    function g = egrad2rgrad(x, g)
        for i = 1 : nelems
            g{i} = elements{i}.egrad2rgrad(...
                                               x{i}, g{i});
        end
    end

    M.ehess2rhess = @ehess2rhess;
    function h = ehess2rhess(x, eg, eh, h)
        for i = 1 : nelems
            h{i} = elements{i}.ehess2rhess(...
                 x{i}, eg{i}, eh{i}, h{i});
        end
    end
    
    M.exp = @exp;
    function y = exp(x, u, t)
        if nargin < 3
            t = 1.0;
        end
        y = cell(nelems, 1);
        for i = 1 : nelems
            y{i} = elements{i}.exp(x{i}, ...
                                                   u{i}, t);
        end
    end
    
    M.retr = @retr;
    function y = retr(x, u, t)
        if nargin < 3
            t = 1.0;
        end
        y = cell(nelems, 1);
        for i = 1 : nelems
            y{i} = elements{i}.retr(x{i}, ...
                                                    u{i}, t);
        end
    end
    
    M.log = @log;
    function u = log(x1, x2)
        u = cell(nelems, 1);
        for i = 1 : nelems
            u{i} = elements{i}.log(x1{i}, ...
                                                   x2{i});
        end
    end

    M.hash = @hash;
    function str = hash(x)
        str = '';
        for i = 1 : nelems
            str = [str elements{i}.hash(x{i})]; %#ok<AGROW>
        end
        str = ['z' hashmd5(str)];
    end

    M.lincomb = @lincomb;
    function v = lincomb(x, a1, u1, a2, u2)
        if nargin == 3
            v = cell(nelems, 1);
            for i = 1 : nelems
                v{i} = elements{i}.lincomb(x{i}, ...
                                                        a1, u1{i});
            end
        elseif nargin == 5
            v = cell(nelems, 1);
            for i = 1 : nelems
                v{i} = elements{i}.lincomb(x{i}, ...
                                     a1, u1{i}, a2, u2{i});
            end
        else
            error('Bad usage of productmanifold.lincomb');
        end
    end

    M.rand = @rand;
    function x = rand()
        x = cell(nelems, 1);
        for i = 1 : nelems
            x{i} = elements{i}.rand();
        end
    end

    M.randvec = @randvec;
    function u = randvec(x)
        u = cell(nelems, 1);
        for i = 1 : nelems
            u{i} = elements{i}.randvec(x{i});
        end
        u = M.lincomb(x, 1/sqrt(nelems), u);
    end

    M.zerovec = @zerovec;
    function u = zerovec(x)
        u = cell(nelems, 1);
        for i = 1 : nelems
            u{i} = elements{i}.zerovec(x{i});
        end
    end

    M.transp = @transp;
    function v = transp(x1, x2, u)
        v = cell(nelems, 1);
        for i = 1 : nelems
            v{i} = elements{i}.transp(x1{i}, ...
                                              x2{i}, u{i});
        end
    end

    M.pairmean = @pairmean;
    function y = pairmean(x1, x2)
        y = cell(nelems, 1);
        for i = 1 : nelems
            y{i} = elements{i}.pairmean(x1{i}, ...
                                                        x2{i});
        end
    end


    % Gather the length of the column vector representations of tangent
    % vectors for each of the manifolds. Raise a flag if any of the base
    % manifolds has no vec function available.
    vec_available = true;
    vec_lens = zeros(nelems, 1);
    for ii = 1 : nelems
        if isfield(elements{ii}, 'vec')
            rand_x = elements{ii}.rand();
            zero_u = elements{ii}.zerovec(rand_x);
            vec_lens(ii) = length(elements{ii}.vec(rand_x, zero_u));
        else
            vec_available = false;
            break;
        end
    end
    vec_pos = cumsum([1 ; vec_lens]);
    
    if vec_available
        M.vec = @vec;
        M.mat = @mat;
    end
    
    function u_vec = vec(x, u_mat)
        u_vec = zeros(vec_pos(end)-1, 1);
        for i = 1 : nelems
            range = vec_pos(i) : (vec_pos(i+1)-1);
            u_vec(range) = elements{i}.vec(x{i}, ...
                                                   u_mat{i});
        end
    end

    function u_mat = mat(x, u_vec)
        u_mat = struct();
        for i = 1 : nelems
            range = vec_pos(i) : (vec_pos(i+1)-1);
            u_mat{i} = elements{i}.mat(x{i}, ...
                                                       u_vec(range));
        end
    end

    vecmatareisometries = true;
    for ii = 1 : nelems
        if ~isfield(elements{ii}, 'vecmatareisometries') || ...
           ~elements{ii}.vecmatareisometries()
            vecmatareisometries = false;
            break;
        end
    end
    M.vecmatareisometries = @() vecmatareisometries;    

end
