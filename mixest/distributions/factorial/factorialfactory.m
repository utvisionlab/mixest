%% |factorialfactory|
% Construct a factorial distribution structure
%
% *Syntax*
%
%   D = factorialfactory(ComponentD, num)
%   D = factorialfactory(ComponentD)
%
% *Description*
%
% |D = factorialfactory(ComponentD, num)| returns a structure representing a
% factorial distribution. |ComponentD| is a distribution structure defining
% facorial distribution type, and |num| is the number of factors.
%
% |D = factorialfactory(ComponentD)| where ComponentD is a cell array of
% distribution structures defined on the same data space, constructs a
% heterogeneous factorial distribution |D| where each component may be of a
% different distribution type.
%
% *Distribution Parameters*
%
% * *|D|* (|num-by-1| cell array of distribution parameter structures) :
% Contains the parameters for each component.
% * *|W|* (|num-by-num| matrix) : The mixing matrix.
%
% *Probability Density Function*
%
% The distribution has the following density:
% 
% $$ f(x)=det(W) \prod_{k=1}^{num} f_k( w_k x) $$
%
% where $num$ is the number of components, $w_k$ is the k-th row of mixing
% matrix $W$, and $f_k$ represents the density function for k-th component.
% 
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

function  D = factorialfactory(ComponentD, num)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'factorial';

%% 
% Flag to control the memory usage (resulting code will be slower)

    low_memory = true;
    
%%

    if nargin > 1
        homogeneous = true;
        Components = num2cell(repmat(ComponentD, [num, 1]));
    else
        homogeneous = false;
        num = numel(ComponentD);
        Components = ComponentD(:);
    end

%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.

    if homogeneous
        ComponentsM = powermanifold(ComponentD.M, num);
    else
        elements = cell(num, 1);
        for k = 1:num %#ok<FXUP>
            elements{k} = Components{k}.M;
        end
        ComponentsM = mxe_productmanifold(elements);
    end
    W = obliquefactory(num, num, true);
    % W = euclideanfactory(num, num);
    D.M = productmanifold(struct('D', ComponentsM, 'W', W));
    
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

%% |component|
% Component distributions
%
% *Syntax*
%
%   D_k = D.component(k)
%

    D.component = @(k) Components{k};

%% |subparam|
% Extract a subset of component parameters
%
% *Syntax*
%
%   subtheta = D.subparam(theta, idx)
%

    D.subparam = @subparam;
    function theta = subparam(theta1, idx)
    % get a subset of parameters from theta1, that are indicated by
    % component indices in idx.
    
        theta.D = theta1.D(idx);
        theta.W = theta1.W(idx,:);
    end

%% |dim|
% See <doc_distribution_common.html#3 distribution structure common members>.

    D.dim = @dim; % parameter space dimensions
    function dim = dim()
        if homogeneous
            dim = ComponentD.dim()*num + num*num;
        else
            dim = num*num; % component weights
            for k = 1:num %#ok<FXUP>
                dim = dim + Components{k}.dim(); % component parameters
            end
        end
    end

%% |datadim|
% See <doc_distribution_common.html#4 distribution structure common members>.

    D.datadim = @() Components{1}.datadim(); % data space dimensions

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
        
        %if ~isfield(store, 'data') || ~isfield(store, 'weight')
        data = mxe_readdata(data);
        data = data.data;
        %end
        
        if ~isfield(store, 'dataTW')
            store.dataTW = theta.W * data;
        end
    end        
    function store = ll_intermediate_params3(weight, store)
    % calculate intermediate parameters for ll #3
        
        if ~isfield(store, 'sumW')
            if ~isempty(weight)
                store.sumW = sum(weight);
            else
                store.sumW = size(store.dataTW,2);
            end
        end
    end

%% |ll|
% See <doc_distribution_common.html#5 distribution structure common members>.

    D.ll = @ll;
    function [ll, store] = ll(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        store = ll_intermediate_params1(theta, store);
        logdetW = store.logdetW;

        store = ll_intermediate_params2(theta, data, store);
        weight = mxe_readweight(data);
        %weight = store.weight;
        dataTW = store.dataTW;
        
        store = ll_intermediate_params3(weight, store);
        sumW = store.sumW;
        
        ll = logdetW * sumW;
        
        if ~isfield(store, 'componentStores')
            store.componentStores = cell(num,1);
        end
        
        for k = 1:num %#ok<FXUP>
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            if low_memory
                lltemp = Components{k}.ll(theta.D{k}, ...
                    struct('data',dataTW(k,:), 'weight',weight));
            else
                [lltemp, store.componentStores{k}] = Components{k}.ll(theta.D{k}, ...
                    struct('data',dataTW(k,:), 'weight',weight), store.componentStores{k});
            end
            ll = ll + lltemp;
        end
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
        
        store = ll_intermediate_params2(theta, data, store);
        weight = mxe_readweight(data);
        %weight = store.weight;
        dataTW = store.dataTW;
        %logdetW = store.logdetW;
        
        llvec = logdetW;
        
        if ~isfield(store, 'componentStores')
            store.componentStores = cell(num,1);
        end
        
        for k = 1:num %#ok<FXUP>
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            if low_memory
                lltemp = Components{k}.llvec(theta.D{k}, ...
                struct('data',dataTW(k,:), 'weight',weight));                
            else
                [lltemp, store.componentStores{k}] = Components{k}.llvec(theta.D{k}, ...
                struct('data',dataTW(k,:), 'weight',weight), store.componentStores{k});
            end
            llvec = llvec + lltemp;
        end
        
        if low_memory
            store = rmfield(store,'dataTW');
            %store = rmfield(store,'componentStores');
        end
    end

%% |llgrad|
% See <doc_distribution_common.html#7 distribution structure common members>.

    D.llgrad = @llgrad;
    function [dll, store] = llgrad(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        store = ll_intermediate_params2(theta, data, store);
        %weight = store.weight;
        dataTW = store.dataTW;
        
        data = mxe_readdata(data);
        weight = data.weight;
        data = data.data;

        store = ll_intermediate_params3(weight, store);
        sumW = store.sumW;
        
        if ~isfield(store, 'componentStores')
            store.componentStores = cell(num,1);
        end
        
        %flag = false;
        if ~isfield(store,'dldD')
            flag = true;
            store.dldD = zeros(size(dataTW));
        end
        for k = 1:num %#ok<FXUP>
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            if flag
                if low_memory
                    dldD = Components{k}.llgraddata(theta.D{k}, ...
                        struct('data',dataTW(k,:), 'weight',weight)); 
                else
                    [dldD, store.componentStores{k}] = ...
                        Components{k}.llgraddata(theta.D{k}, ...
                        struct('data',dataTW(k,:), 'weight',weight), ...
                        store.componentStores{k});
                end
                store.dldD(k,:) = dldD;
            end
            if low_memory
                dll.D{k} = Components{k}.llgrad(theta.D{k}, ...
                    struct('data',dataTW(k,:), 'weight',weight));
            else
                [dll.D{k}, store.componentStores{k}] = ...
                    Components{k}.llgrad(theta.D{k}, ...
                    struct('data',dataTW(k,:), 'weight',weight), ...
                    store.componentStores{k});
            end
            % dll.W(k,:) = sum(bsxfun(@times, dldD, data), 2).';
        end   
        dll.W = store.dldD * data.';
        dll.W = dll.W + (sumW) * inv(theta.W).';
        
        if low_memory
            store = rmfield(store,'dataTW');
            store = rmfield(store,'dldD');
            %store = rmfield(store,'componentStores');
        end
    end

%% |llgraddata|
% See <doc_distribution_common.html#8 distribution structure common members>.

    D.llgraddata = @llgraddata;
    function [dld, store] = llgraddata(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        store = ll_intermediate_params2(theta, data, store);
        weight = mxe_readweight(data);
        %weight = store.weight;
        dataTW = store.dataTW;

        if ~isfield(store, 'componentStores')
            store.componentStores = cell(num,1);
        end
        
        dld = 0;
        store.dldD = zeros(size(dataTW));

        for k = 1:num %#ok<FXUP>
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            if low_memory
                dldD = Components{k}.llgraddata(theta.D{k}, ...
                    struct('data',dataTW(k,:), 'weight',weight));
            else
                [dldD, store.componentStores{k}] = ...
                    Components{k}.llgraddata(theta.D{k}, ...
                    struct('data',dataTW(k,:), 'weight',weight), ...
                    store.componentStores{k});
            end
            store.dldD(k,:) = dldD;
            dld = dld + bsxfun(@times, dldD, theta.W(k,:).');  
        end  
        
        if low_memory
            store = rmfield(store,'dataTW');
            store = rmfield(store,'dldD');
            %store = rmfield(store,'componentStores');
        end
    end

%% |gaussianize|
% 
    D.gaussianize = @gaussianize;
    function y = gaussianize(theta, data)
        data = mxe_readdata(data);
        y = theta.W * data.data;
        for k = 1:num %#ok<FXUP>
            y (k, :) = Components{k}.gaussianize(theta.D{k}, y(k,:));
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
        for k = 1:num %#ok<FXUP>
            data (k, :) = Components{k}.sample(theta.D{k}, n);
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
        for k = 1:num %#ok<FXUP>
            theta.D{k} = Components{k}.init(data(k,:), varargin{:});
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
        data = data.data(:).';
        for k = 1:num %#ok<FXUP>
            penalizer_theta.D{k} = Components{k}.penalizerparam(data);
        end
    end

%% |penalizercost|
% See <doc_distribution_common.html#16 distribution structure common members>.

    D.penalizercost = @penalizercost;
    function [costP, store] = penalizercost(theta, penalizer_theta, store)
        
        if nargin < 3
            store = struct;
        end
        
        if ~isfield(store, 'componentStores')
            store.componentStores = cell(num,1);
        end
        
        costP = 0;
        for k = 1:num %#ok<FXUP>
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            [cost, store.componentStores{k}] = ...
                Components{k}.penalizercost(theta.D{k}, penalizer_theta.D{k}, ...
                store.componentStores{k});
            costP = costP + cost;
        end
    end

%% |penalizergrad|
% See <doc_distribution_common.html#17 distribution structure common members>.

    D.penalizergrad = @penalizergrad;
    function [gradP, store] = penalizergrad(theta, penalizer_theta, store)
        
        if nargin < 3
            store = struct;
        end
        
        if ~isfield(store, 'componentStores')
            store.componentStores = cell(num,1);
        end
        
        gradP.W = zeros(num,num);
        for k = 1:num %#ok<FXUP>
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            [gradP.D{k}, store.componentStores{k}] = ...
                Components{k}.penalizergrad(theta.D{k}, penalizer_theta.D{k}, ...
                store.componentStores{k});
        end
    end

%% |sumparam|
% See <doc_distribution_common.html#18 distribution structure common members>.

    D.sumparam = @sumparam;
    function theta = sumparam(theta1, theta2)
        theta.D = cell(num, 1);
        for k = 1:num %#ok<FXUP>
            theta.D{k} = Components{k}.sumparam(theta1.D{k}, theta2.D{k});
        end
        theta.W = theta1.W + theta2.W;
    end

%% |scaleparam|
% See <doc_distribution_common.html#19 distribution structure common members>.

    D.scaleparam = @scaleparam;
    function theta = scaleparam(scalar, theta)
        for k = 1:num %#ok<FXUP>
            theta.D{k} = Components{k}.scaleparam(scalar, theta.D{k});
        end
        theta.W = scalar * theta.W;
    end

%% |sumgrad|
% See <doc_distribution_common.html#20 distribution structure common members>.

%% |scalegrad|
% See <doc_distribution_common.html#21 distribution structure common members>.

%% |entropy|
% See <doc_distribution_common.html#22 distribution structure common members>.

    D.entropy = @entropy;
    function h = entropy(theta)
        h = zeros(1,num);
        for k = 1:num %#ok<FXUP>
            h(k) = Components{k}.entropy(theta.D{k});
        end
    end

%% |kl|
% See <doc_distribution_common.html#23 distribution structure common members>.

    D.kl = @kl;
    function kl = kl(theta)
        assert(homogeneous, 'KL-divergence is only valid for homogeneous mixtures.');
        kl = zeros(num);
        for k1 = 1:num
            for k2 = 1:k1-1
                kl(k1,k2) = ComponentD.kl(theta.D{k1}, theta.D{k2});
                kl(k2,k1) = kl(k1,k2);
            end
        end        
    end
    
%% |AICc|
% See <doc_distribution_common.html#24 distribution structure common members>.

%% |BIC|
% See <doc_distribution_common.html#25 distribution structure common members>.

%% |display|
% See <doc_distribution_common.html#26 distribution structure common members>.

    D.display = @display;
    function str = display(theta)
        str = '';
        for k = 1:num %#ok<FXUP>
            str = [str, sprintf('\ncomponent(%d): %s / w_k: %s\n', k, Components{k}.name(), mat2str(theta.W(k,:), 4))]; %#ok<AGROW>
            str = [str, Components{k}.display(theta.D{k})]; %#ok<AGROW>
        end
        
        if nargout == 0
            str = [sprintf('%s distribution parameters:\n', D.name()), str];
        end
    end

%%

    D = mxe_addsharedfields(D);
end
