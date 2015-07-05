%% |factorialmogfactory|
% Construct a factorial distribution with MoG distributions for factors
%
% *Syntax*
%
%   D = factorialmogfactory(numCpt, num)
%
% *Description*
%
% |D = factorialfactory(numCpt, num)| returns a structure representing a
% factorial distribution. |numCpt| is number of Gaussian components
% for factors, and |num| is the number of factors.
%
% *Distribution Parameters*
%
% * *|mu|* (|numCpt-by-num| matrix) : Mean values
% * *|sd|* (|numCpt-by-num| matrix) : Standard deviations
% * *|W|* (|num-by-num| matrix) : The mixing matrix.
%
% *Probability Density Function*
%
% The distribution has the following density:
% 
% $$ f(x)=det(W) \prod_{k=1}^{num} f_k( w_k x) $$
%
% where $num$ is the number of components, $w_k$ is the k-th row of mixing
% matrix $W$, and $f_k$ is a MoG the density function for k-th component.
% 
% Important: When factors are MoG, it is faster to use this factory instead 
% of factorialfactory 
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Contributors:
%  Reshad Hosseini
%  Mohamadreza Mash'al
%
% Change log: 
%

function  D = factorialmogfactory(numCpt, num)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'factorial';

%% 
% Flag to control the memory usage (resulting code will be slower)

    low_memory = true;

%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.

    W = obliquefactory(num, num, true);
    mu = euclideanfactory(numCpt, num);
    sd = euclideanfactory(numCpt, num);
    p = productsimplexfactory(numCpt, num);
    D.M = productmanifold(struct('p', p, 'mu', mu, 'sd', sd, 'W', W));
    
    % from here on use Components{k} to access k'th component irrespective
    % to the mixture being homogeneous or heterogeneous.

%% |num|
% Number of components (excluding any fixed components)
%
% *Syntax*
%
%   num = D.num()
%

    D.num = @() num;

%% |numCpt|
% Component distributions
%
% *Syntax*
%
%   numCpt = D.numCpt()
%

    D.numCpt = @() numCpt;

%% |dim|
% See <doc_distribution_common.html#3 distribution structure common members>.

    D.dim = @dim; % parameter space dimensions
    function dim = dim()
        dim = (3*numCpt-1)*num + (num-1)*num;
    end

%% |datadim|
% See <doc_distribution_common.html#4 distribution structure common members>.

    D.datadim = @() num(); % data space dimensions

%%

    function store = ll_intermediate_params1(theta, store)
    % calculate intermediate parameters for ll #2
        
        if ~isfield(store, 'logdetW')
            % LU factorization needed for computing Log det
            Y = lu(theta.W);
            store.logdetW = sum(log(abs(diag(Y))));
            %logdetSign = prod(sign(diag(Y)))
        end
    end
    function store = ll_intermediate_params2(theta, data, store)
    % calculate intermediate parameters for ll #1
        
        if ~isfield(store, 'dataTW')
            data = mxe_readdata(data);
            data = data.data;
            store.dataTW = theta.W * data;
        end
    end        
    function store = weighting_intermediate_params(theta, data, store)
    % calculate intermediate parameters for ll #4
        
        if ~isfield(store, 'llse') || ~isfield(store, 'llvec')
            store = ll_intermediate_params2(theta, data, store);
            % computing log-likelihood of all components
            lls = repmat(store.dataTW, [ 1 1 numCpt]);

            lls = bsxfun(@minus, lls, shiftdim(shiftdim(theta.mu, -1),2));
            lls = bsxfun(@rdivide, lls, shiftdim(shiftdim(theta.sd, -1),2));
            
            cte = -0.5*log(2*pi) - 0.5*shiftdim(shiftdim(log(theta.sd.^2), -1),2) + ...
                shiftdim(shiftdim(log(theta.p), -1),2);
            lls = bsxfun(@plus, -0.5 * lls.^2, cte);
            
            % using lls compute weights for mixtures
            store.llvec = logsumexp(lls, 3);
            store.lls = lls;
        end
    end
    function [component_weights, store] = weighting(theta, data, store)
    % calculate weighting
        
        if nargin < 3
            store = struct;
        end
        
        store = weighting_intermediate_params(theta, data, store);
        lls = store.lls;
        llvec = store.llvec;
        
        weight = mxe_readweight(data);
            
        % undo logarithm and normalize hX such that each column sums up to 1.
        component_weights = exp( bsxfun(@minus, lls, llvec) );
        
        if ~isempty(weight)
            component_weights = bsxfun(@times, component_weights, weight);
        end
    end

%% |ll|
% See <doc_distribution_common.html#5 distribution structure common members>.

    D.ll = @ll;
    function [ll, store] = ll(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        [llvecc, store] = llvec(theta, data, store);
        
        ll = sum (llvecc, 2);

    end

%% |llvec|
% See <doc_distribution_common.html#6 distribution structure common members>.

    D.llvec = @llvec;
    function [llvec, store] = llvec(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        store = ll_intermediate_params1(theta, store);
        logdetW = store.logdetW;
        
        store = weighting_intermediate_params(theta, data, store);
        weight = mxe_readweight(data);
        
        llvec = logdetW;

        llvec = llvec + sum(store.llvec);
        
        if ~isempty(weight)
            llvec = llvec .* weight;
        end
        
        if low_memory
            store = rmfield(store,'dataTW');
            store = rmfield(store,'lls');
            store = rmfield(store,'llvec');
        end
    end

%% |llgrad|
% See <doc_distribution_common.html#7 distribution structure common members>.

    D.llgrad = @llgrad;
    function [dll, store] = llgrad(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        [component_weights, store] = weighting(theta, data, store);
        dataTW = store.dataTW;
        
        data = mxe_readdata(data);
        weight = data.weight;
        data = data.data;
        
        dldD = repmat(dataTW, [ 1 1 numCpt]);
        dldD = bsxfun(@minus, dldD, shiftdim(shiftdim(theta.mu, -1),2));
        dldD = bsxfun(@rdivide, dldD, shiftdim(shiftdim(theta.sd, -1),2).^2);
        dldD = component_weights.*dldD;
        dll.mu = squeeze(sum(dldD,2)).';
        store.dldD = squeeze(sum(-dldD,3));
        dldD = repmat(dataTW, [ 1 1 numCpt]);
        dldD = bsxfun(@minus, dldD, shiftdim(shiftdim(theta.mu, -1),2));
        dldD = bsxfun(@rdivide, dldD.^2, shiftdim(shiftdim(theta.sd.^3, -1),2));
        dldD = component_weights.*dldD;
        dll.sd = squeeze(sum(dldD,2));
        sumW = squeeze(sum(component_weights, 2)).';
        dll.sd = dll.sd.' - sumW./theta.sd;
        dll.p = sumW ./ theta.p;
 
        if isempty(weight)
            weight = size(store.dldD, 2);
        end
        
        dll.W = store.dldD * data.';
        dll.W = dll.W + sum(weight) * inv(theta.W).';
        
        if low_memory
            store = rmfield(store,'dataTW');
            store = rmfield(store,'dldD');
            store = rmfield(store,'lls');
            store = rmfield(store,'llvec');
        end
    end

%% |llgraddata|
% See <doc_distribution_common.html#8 distribution structure common members>.

    D.llgraddata = @llgraddata;
    function [dld, store] = llgraddata(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        if nargin < 3
            store = struct;
        end
        
        [component_weights, store] = weighting(theta, data, store);
        dataTW = store.dataTW;
        
        if ~isfield('store','dldD')
            dldD = repmat(dataTW, [ 1 1 numCpt]);
            dldD = bsxfun(@minus, dldD, shiftdim(shiftdim(theta.mu, -1),2));
            dldD = bsxfun(@rdivide, dldD, shiftdim(shiftdim(theta.sd, -1),2).^2);
            dldD = component_weights.*dldD;
            store.dldD = squeeze(sum(-dldD,3));
        end
        
        dld = 0;
        for k = 1:num
            dld = dld + bsxfun(@times, store.dldD(k,:), theta.W(k,:).'); 
        end
        
        if low_memory
            store = rmfield(store,'dataTW');
            store = rmfield(store,'dldD');
            store = rmfield(store,'lls');
            store = rmfield(store,'llvec');
        end
    end

%% |gaussianize|
% 
    D.gaussianize = @gaussianize;
    function y = gaussianize(theta, data)
        data = mxe_readdata(data);
        y = theta.W * data.data;
        Components = mixturefactory(mvnfactory(1),numCpt);
        for k = 1:num %#ok<FXUP>
            for l = 1:numCpt
                t.D{l}.mu = theta.mu(l,k);
                t.D{l}.sigma = theta.sd(l,k).^2;
            end
            t.p = theta.p(:,k);
            y (k, :) = Components.gaussianize(t, y(k,:));
        end
    end

%% |pdf|
% See <doc_distribution_common.html#10 distribution structure common members>.

%% |sample|
% See <doc_distribution_common.html#11 distribution structure common members>.

    D.sample = @sample;
    function data = sample(theta, n)
        
        if nargin < 2, n = 1; end
        
        data = zeros(num, n);
        Components = mixturefactory(mvnfactory(1),numCpt);
        for k = 1:num %#ok<FXUP>
            for l = 1:numCpt
                t.D{l}.mu = theta.mu(l,k);
                t.D{l}.sigma = theta.sd(l,k).^2;
            end
            t.p = theta.p(:,k);
            data (k, :) = Components.sample(t, n);
        end
        data = theta.W \ data;
    end

%% |randparam|
% See <doc_distribution_common.html#12 distribution structure common members>.

%% |init|
% See <doc_distribution_common.html#13 distribution structure common members>.

    D.init = @init;
    function theta = init(data, varargin)
        data = mxe_readdata(data);
        data = data.data;
        
        theta = D.M.rand();
        data = theta.W * data;
        Components = mixturefactory(mvnfactory(1),numCpt);
        for k = 1:num %#ok<FXUP>
            t = Components.init(data(k,:), varargin{:});
            for l = 1:numCpt
                theta.mu(l,k) = t.D{l}.mu;
                theta.sd(l,k) = sqrt(t.D{l}.sigma);
            end
            theta.p(:,k) = t.p;
        end
    end

%%

% %% |estimatedefault|
%     D.estimatedefault = @estimatedefault;
%     function [varargout] = estimatedefault(varargin)
%         [varargout{1:nargout}] = factorial_estimatedefault(D, varargin{:}); 
%     end

%% |penalizerparam|
% See <doc_distribution_common.html#15 distribution structure common members>.
%
% *Penalizer Info*
%
% The default penalizer for the mixture distribution is the mixture of the
% default penalizers of its components, with equal weights.
%

%TODO doc

    D.penalizerparam = @penalizerparam;
    function penalizer_theta = penalizerparam(data)
        data = mxe_readdata(data);
        n = data.size;
        data = data.data;

        penalizer_theta.kappa = 0.01; %0.1;
        penalizer_theta.nu = 3; %dim+2; %2;
        penalizer_theta.mu = sum(data(:))/n/num;
        data = bsxfun(@minus, data(:), penalizer_theta.mu);
        sigma = (data.' * data)/n/num;
        penalizer_theta.invLambda = sigma / numCpt^2;
        
    end

%% |penalizercost|
% See <doc_distribution_common.html#16 distribution structure common members>.

    D.penalizercost = @penalizercost;
    function [costP, store] = penalizercost(theta, penalizer_theta, store)
        
        if nargin < 3
            store = struct;
        end
        
        logdetR = log(theta.sd);
        Sinv = 1./ theta.sd.^2;
        const1 = 0.5 * penalizer_theta.kappa;
        const2 = penalizer_theta.nu + 2;
        
        mu = theta.mu - penalizer_theta.mu;
        SinvMu = Sinv .* mu;
        
        if penalizer_theta.kappa == 0
            costP = 0;
        else
            costP = - logdetR - const1 * mu .* SinvMu;
        end
        
        % If data to penalizeparam is whitened then invLambda is scalar
        costP = costP - const2 * logdetR - 0.5 * penalizer_theta.invLambda * Sinv;
        
        costP = sum(costP(:));
        
    end

%% |penalizergrad|
% See <doc_distribution_common.html#17 distribution structure common members>.

    D.penalizergrad = @penalizergrad;
    function [gradP, store] = penalizergrad(theta, penalizer_theta, store)
        
        if nargin < 3
            store = struct;
        end
        
        Sinv = 1./ theta.sd.^2;
        const1 = 0.5 * penalizer_theta.kappa;
        const2 = penalizer_theta.nu + 2;
        
        if penalizer_theta.kappa == 0
            const2 = - 0.5 * const2;
        else
            const2 = - 0.5 * (const2 + 1);
        end
        
        mu = theta.mu - penalizer_theta.mu;
        SinvMu = Sinv .* mu;
        
        
        
        gradP.sd = const2 * Sinv + const1 * (SinvMu .* SinvMu) + ...
            (0.5 * penalizer_theta.invLambda) * (Sinv .* Sinv);
        
        gradP.sd = 2 * gradP.sd .* theta.sd;

        gradP.mu = - penalizer_theta.kappa * SinvMu;
        
        gradP.p = zeros(size(theta.p));
        gradP.W = zeros(size(theta.W));
        
    end

%% |sumparam|
% See <doc_distribution_common.html#18 distribution structure common members>.

    D.sumparam = @sumparam;
    function theta = sumparam(theta1, theta2)
        theta.p = theta1.p + theta2.p;
        theta.mu = theta1.mu + theta2.mu;
        theta.sd = theta1.sd + theta2.sd;
        theta.W = theta1.W + theta2.W;
    end

%% |scaleparam|
% See <doc_distribution_common.html#19 distribution structure common members>.

    D.scaleparam = @scaleparam;
    function theta = scaleparam(scalar, theta)
        theta.p = scalar * theta.p;
        theta.mu = scalar * theta.mu;
        theta.sd = scalar * theta.sd;
        theta.W = scalar * theta.W;
    end

%% |sumgrad|
% See <doc_distribution_common.html#20 distribution structure common members>.

%% |scalegrad|
% See <doc_distribution_common.html#21 distribution structure common members>.

%% |entropy|
% See <doc_distribution_common.html#22 distribution structure common members>.

%% |kl|
% See <doc_distribution_common.html#23 distribution structure common members>.

    
%% |AICc|
% See <doc_distribution_common.html#24 distribution structure common members>.

%% |BIC|
% See <doc_distribution_common.html#25 distribution structure common members>.

%% |display|
% See <doc_distribution_common.html#26 distribution structure common members>.

    D.display = @display;
    function str = display(theta)
        str = '';
         Components = mixturefactory(mvnfactory(1),numCpt);        
        for k = 1:num %#ok<FXUP>
            str = [str, sprintf('\ncomponent(%d): %s / w_k: %s\n', k, Components.name(), mat2str(theta.W(k,:), 4))]; %#ok<AGROW>
            for l = 1:numCpt
                t.D{l}.mu = theta.mu(l,k);
                t.D{l}.sigma = theta.sd(l,k).^2;
            end
            t.p = theta.p(:,k);
            str = [str, Components.display(t)]; %#ok<AGROW>
        end
        
        if nargout == 0
            str = [sprintf('%s distribution parameters:\n', D.name()), str];
        end
    end

%%

    D = mxe_addsharedfields(D);
end
