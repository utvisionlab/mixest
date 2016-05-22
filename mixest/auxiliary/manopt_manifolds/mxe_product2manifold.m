function M = mxe_product2manifold(elements)
% 
%
        
        
        
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
            val = val + elements{i}.inner(x.(elems{i}), ...
                                               u.(elems{i}), v.(elems{i}));
        end
    end

    M.norm = @(x, d) sqrt(M.inner(x, d, d));

    M.dist = @dist;
    function d = dist(x, y)
        sqd = 0;
        for i = 1 : nelems
            sqd = sqd + elements{i}.dist(x.(elems{i}), ...
                                                 y.(elems{i}))^2;
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
        for i = 1 : nelems
            v.(elems{i}) = elements{i}.proj(x.(elems{i}), ...
                                                    u.(elems{i}));
        end
    end

    M.tangent = @tangent;
    function v = tangent(x, u)
        for i = 1 : nelems
            v.(elems{i}) = elements{i}.tangent(x.(elems{i}), ...
                                                       u.(elems{i}));
        end
    end

    M.tangent2ambient = @tangent2ambient;
    function v = tangent2ambient(x, u)
        for i = 1 : nelems
            if isfield(elements{i}, 'tangent2ambient')
                v.(elems{i}) = ...
                    elements{i}.tangent2ambient( ...
                                               x.(elems{i}), u.(elems{i}));
            else
                v.(elems{i}) = u.(elems{i});
            end
        end
    end

    M.egrad2rgrad = @egrad2rgrad;
    function g = egrad2rgrad(x, g)
        for i = 1 : nelems
            g.(elems{i}) = elements{i}.egrad2rgrad(...
                                               x.(elems{i}), g.(elems{i}));
        end
    end

    M.ehess2rhess = @ehess2rhess;
    function h = ehess2rhess(x, eg, eh, h)
        for i = 1 : nelems
            h.(elems{i}) = elements{i}.ehess2rhess(...
                 x.(elems{i}), eg.(elems{i}), eh.(elems{i}), h.(elems{i}));
        end
    end
    
    M.exp = @exp;
    function y = exp(x, u, t)
        if nargin < 3
            t = 1.0;
        end
        for i = 1 : nelems
            y.(elems{i}) = elements{i}.exp(x.(elems{i}), ...
                                                   u.(elems{i}), t);
        end
    end
    
    M.retr = @retr;
    function [y, store] = retr(x, u, t, store)
        if nargin < 3
            t = 1.0;
        end
        if nargin < 4 && nargout == 1
            for i = 1 : nelems
                y.(elems{i}) = elements{i}.retr(x.(elems{i}), ...
                    u.(elems{i}), t);
            end
        end
        if nargin == 4 && nargout == 1
            for i = 1 : nelems
                y.(elems{i}) = elements{i}.retr(x.(elems{i}), ...
                    u.(elems{i}), t, store.(elems{i}));
            end
        end
        if nargin < 4 && nargout == 2
            for i = 1 : nelems
                [y.(elems{i}), store.(elems{i})] = elements{i}.retr(x.(elems{i}), ...
                    u.(elems{i}), t);
            end
        end
        if nargin == 4 && nargout == 2
            for i = 1 : nelems
                [y.(elems{i}), store.(elems{i})] = elements{i}.retr(x.(elems{i}), ...
                    u.(elems{i}), t, store.(elems{i}));
            end
        end
    end
    
    M.log = @log;
    function u = log(x1, x2)
        for i = 1 : nelems
            u.(elems{i}) = elements{i}.log(x1.(elems{i}), ...
                                                   x2.(elems{i}));
        end
    end

    M.hash = @hash;
    function str = hash(x)
        str = '';
        for i = 1 : nelems
            str = [str elements{i}.hash(x.(elems{i}))]; %#ok<AGROW>
        end
        str = ['z' hashmd5(str)];
    end

    M.lincomb = @lincomb;
    function v = lincomb(x, a1, u1, a2, u2)
        if nargin == 3
            for i = 1 : nelems
                v.(elems{i}) = elements{i}.lincomb(x.(elems{i}), ...
                                                        a1, u1.(elems{i}));
            end
        elseif nargin == 5
            for i = 1 : nelems
                v.(elems{i}) = elements{i}.lincomb(x.(elems{i}), ...
                                     a1, u1.(elems{i}), a2, u2.(elems{i}));
            end
        else
            error('Bad usage of productmanifold.lincomb');
        end
    end

    M.rand = @rand;
    function x = rand()
        for i = 1 : nelems
            x.(elems{i}) = elements{i}.rand();
        end
    end

    M.randvec = @randvec;
    function u = randvec(x)
        for i = 1 : nelems
            u.(elems{i}) = elements{i}.randvec(x.(elems{i}));
        end
        u = M.lincomb(x, 1/sqrt(nelems), u);
    end

    M.zerovec = @zerovec;
    function u = zerovec(x)
        for i = 1 : nelems
            u.(elems{i}) = elements{i}.zerovec(x.(elems{i}));
        end
    end

    M.transp = @transp;
    function [v, store] = transp(x1, x2, e, u, t, store)
        if nargin == 5 && nargout == 1
            for i = 1 : nelems
                v.(elems{i}) = elements{i}.transp(x1.(elems{i}), ...
                    x2.(elems{i}), e.(elems{i}), u.(elems{i}), t);
            end
        end
        if nargin == 3 && nargout == 1
            for i = 1 : nelems
                v.(elems{i}) = elements{i}.transp(x1.(elems{i}), ...
                    x2.(elems{i}), e.(elems{i}));
            end
        end
        if nargin == 6 && nargout == 1
            for i = 1 : nelems
                v.(elems{i}) = elements{i}.transp(x1.(elems{i}), ...
                    x2.(elems{i}), e.(elems{i}), u.(elems{i}), t, store.(elems{i}));
            end
        end
        if nargin == 5 && nargout == 2
            for i = 1 : nelems
                [v.(elems{i}), store.(elems{i})] = elements{i}.transp(x1.(elems{i}), ...
                    x2.(elems{i}), e.(elems{i}), u.(elems{i}), t);
            end
        end
        if nargin == 3 && nargout == 2
            for i = 1 : nelems
                [v.(elems{i}), store.(elems{i})] = elements{i}.transp(x1.(elems{i}), ...
                    x2.(elems{i}), e.(elems{i}));
            end
        end
        if nargin == 6 && nargout == 2
            for i = 1 : nelems
                [v.(elems{i}), store.(elems{i})] = elements{i}.transp(x1.(elems{i}), ...
                    x2.(elems{i}), e.(elems{i}), u.(elems{i}), t, store.(elems{i}));
            end
        end
    end

    M.retrtransp = @retrtransp;
    function [y, v, store] = retrtransp(x, u, e, t, store)
        if nargin < 3
            t = 1.0;
        end
        if nargin < 5 && nargout == 2
            for i = 1 : nelems
                [y.(elems{i}), v.(elems{i})] = elements{i}.retrtransp(x.(elems{i}), ...
                    u.(elems{i}), e.(elems{i}), t);
            end
        end
        if nargin == 5 && nargout == 2
            for i = 1 : nelems
                [y.(elems{i}), v.(elems{i})] = elements{i}.retrtransp(x.(elems{i}), ...
                    u.(elems{i}), e.(elems{i}), t, store.(elems{i}));
            end
        end
        if nargin < 5 && nargout == 3
            for i = 1 : nelems
                [y.(elems{i}), v.(elems{i}), store.(elems{i})] = elements{i}.retrtransp(x.(elems{i}), ...
                    u.(elems{i}), e.(elems{i}), t);
            end
        end
        if nargin == 5 && nargout == 3
            for i = 1 : nelems
                [y.(elems{i}), v.(elems{i}), store.(elems{i})] = elements{i}.retrtransp(x.(elems{i}), ...
                    u.(elems{i}), e.(elems{i}), t, store.(elems{i}));
            end
        end
    end

    M.transpdiffE = @transpdiffE;
    function [v, store] = transpdiffE(x, u, e, t, store)
        if nargin < 3
            t = 1.0;
        end
        if nargin < 5 && nargout == 1
            for i = 1 : nelems
                [v.(elems{i})] = elements{i}.transpdiffE(x.(elems{i}), ...
                    u.(elems{i}), e.(elems{i}), t);
            end
        end
        if nargin == 5 && nargout == 1
            for i = 1 : nelems
                [v.(elems{i})] = elements{i}.transpdiffE(x.(elems{i}), ...
                    u.(elems{i}), e.(elems{i}), t, store.(elems{i}));
            end
        end
        if nargin < 5 && nargout == 2
            for i = 1 : nelems
                [v.(elems{i}), store.(elems{i})] = elements{i}.transpdiffE(x.(elems{i}), ...
                    u.(elems{i}), e.(elems{i}), t);
            end
        end
        if nargin == 5 && nargout == 2
            for i = 1 : nelems
                [v.(elems{i}), store.(elems{i})] = elements{i}.transpdiffE(x.(elems{i}), ...
                    u.(elems{i}), e.(elems{i}), t, store.(elems{i}));
            end
        end
    end

    M.transpstore = @transpstore;
    function [ec, iec] = transpstore(x, y)
        for i = 1 : nelems
            [ec.(elems{i}), iec.(elems{i})]  = ...
                elements{i}.transpstore(x.(elems{i}), ...
                y.(elems{i}));
        end
    end   
    
    M.transpf = @transpvecfast;
    function y = transpvecfast(x, ec)
        for i = 1 : nelems
            y.(elems{i}) = elements{i}.transpf(x.(elems{i}), ...
                                                        ec.(elems{i}));
        end
    end
    
    M.atranspf = @atranspvecfast;
    function y = atranspvecfast(x, iec)
        for i = 1 : nelems
            y.(elems{i}) = elements{i}.atranspf(x.(elems{i}), ...
                                                        iec.(elems{i}));
        end
    end
    
    M.pairmean = @pairmean;
    function y = pairmean(x1, x2)
        for i = 1 : nelems
            y.(elems{i}) = elements{i}.pairmean(x1.(elems{i}), ...
                                                        x2.(elems{i}));
        end
    end


    % Gather the length of the column vector representations of tangent
    % vectors for each of the manifolds. Raise a flag if any of the base
    % manifolds has no vec function available.
    vec_available = true;
    vec_lens = zeros(nelems, 1);
    for ii = 1 : nelems
        Mi = elements.(elems{ii});
        if isfield(Mi, 'vec')
            rand_x = Mi.rand();
            zero_u = Mi.zerovec(rand_x);
            vec_lens(ii) = length(Mi.vec(rand_x, zero_u));
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
            u_vec(range) = elements{i}.vec(x.(elems{i}), ...
                                                   u_mat.(elems{i}));
        end
    end

    function u_mat = mat(x, u_vec)
        u_mat = struct();
        for i = 1 : nelems
            range = vec_pos(i) : (vec_pos(i+1)-1);
            u_mat.(elems{i}) = elements{i}.mat(x.(elems{i}), ...
                                                       u_vec(range));
        end
    end

    vecmatareisometries = true;
    for ii = 1 : nelems
        if ~isfield(elements.(elems{ii}), 'vecmatareisometries') || ...
           ~elements.(elems{ii}).vecmatareisometries()
            vecmatareisometries = false;
            break;
        end
    end
    M.vecmatareisometries = @() vecmatareisometries;    

end
