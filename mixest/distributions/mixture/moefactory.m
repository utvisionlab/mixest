%% |moefactory|
% Construct a mixture of experts distribution structure
%
% *Syntax*
%
%   D = mixturefactory(GateD, ExpertD)
%
% *Description*
%
% |D = mixturefactory(GateD, ExpertD)| returns a structure representing a
% mixture of experts distribution. |GatetD| is a gate structure. |ExpertD|
% defines the expert distributions and it may be either a conditional
% distribution structure or a cell array of conditional distribution
% structures on the same data space.
%
% When |ExpertD| is a conditional distribution structure, it defines the
% expert component distribution type. The number of experts equals
% |GateD.num()|.
%
% When |ExpertD| is a cell array of conditional distribution structures,
% the function constructs a heterogeneous mixture of experts distribution
% where each expert may be of a different distribution type. The number of
% elements in the cell array should equal |GateD.num()|.
%
% *Distribution Parameters*
%
% * *|D|* (|num-by-1| cell array of distribution parameter structures) :
% Contains the parameters for each expert.
% * *|G|* : Sturucture containing the parameters of the gate.
%
% *Probability Density Function*
%
% The distribution has the following density:
% 
% $$ f(y | x) =  \sum_{k=1}^{num}{ \pi_k(x) g_k(y|x)} $$
%
% where $\pi_k(x)$ is the k-th gate component defined by the gate
% distribution and $g_k(y|x)$ is the k-th expert component.
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

function  D = moefactory(GateD, ExpertD)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'moefactory';

%%
    num = GateD.num();
    datadim = GateD.datadim();

%%
    if ~iscell(ExpertD)
        homogeneous = true;
        Experts = num2cell(repmat(ExpertD, [num,1]));
    else
        homogeneous = false;
        Experts = ExpertD(:);
    end
    
%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.
    if homogeneous
        ExpertM = powermanifold(ExpertD.M, num);
    else
        elements = cell(num, 1);
        for kk = 1:num
            elements{kk} = Experts{kk}.M;
        end
        ExpertM = mxe_productmanifold(elements);
    end
    D.M = productmanifold(struct('D', ExpertM, 'G', GateD.M));

%% |expert|
% expert distributions
%
% *Syntax*
%
%   D_k = D.expert(idx)
%
% *Description*
%
% |D_k = D.expert(idx)| where |idx| is the index of a component, gives
% the component distribution structure at |idx|.
%
% |D_k = D.expert(idx)| where |idx| is an index vector with more than
% one element, returns a cell array of the distribution structures indexed
% in |idx|.
%
% Valid range for the indices is from 1 to |D.numtotal()|.
%

    D.expert = @getexpert;
    function result = getexpert(idx)
        if numel(idx) == 1
            result = Experts{idx};
        else
            result = Experts(idx);
        end
    end

%% |gate|
% Gate distribution

    D.gate = @() GateD;

%% |num|
% Number of components (excluding any fixed components)
%
% *Syntax*
%
%   num = D.num()
%

    D.num = @() num;

%% |dim|
% See <doc_distribution_common.html#3 distribution structure common members>.

    D.dim = @dim; % parameter space dimensions
    function dim = dim()
        dim = GateD.dim();
        for k = 1:num
            dim = dim + Experts{k}.dim(); % component parameters
        end
    end

%% |datadim|
% See <doc_distribution_common.html#4 distribution structure common members>.

    D.datadim = @() datadim(); % data space dimensions

%%

    function store = weighting_intermediate_params(theta, data, store)
    % calculate intermediate parameters for weighting
        
        if ~isfield(store, 'componentStores')
            store.componentStores = cell(num,1);
        end
        
        if ~isfield(store, 'hX') || ~isfield(store, 'llik')
            data = mxe_readdata(data, false);
            index = data.index;
            n = data.size;
            datamat = data.data;
            
            % Initialize different variables
            store.hX = zeros(num, n);
            % Calculate the log-likelihood of data goes to different clusters
            for k = 1:num
                if ~isstruct(store.componentStores{k})
                    store.componentStores{k} = struct;
                end
                [llvec, store.componentStores{k}] = ...
                    Experts{k}.llvec(...
                    theta.D{k}, struct('data',datamat, 'weight',[], 'index',index), ...
                    store.componentStores{k});
                store.hX(k,:) = llvec;
            end
        end
    end

%% |llvec|
% See <doc_distribution_common.html#6 distribution structure common members>.

    D.llvec = @llvec;
    function [llvec, store] = llvec(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        store = weighting_intermediate_params(theta, data, store);

        data = mxe_readdata(data, false);
        
        data.data = [data.data(1:datadim, data.index); store.hX];
        data.index = [];
        [llvec, store] = GateD.llvec(theta.G, data, store);
        
    end

%% |ll|
% See <doc_distribution_common.html#5 distribution structure common members>.

    D.ll = @ll;
    function [ll, store] = ll(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        store = weighting_intermediate_params(theta, data, store);
        
        data = mxe_readdata(data);
        
        data.data = [data.data(1:datadim,:); store.hX];
        
        [ll, store] = GateD.ll(theta.G, data, store);

    end

%% |llgrad|
% See <doc_distribution_common.html#7 distribution structure common members>.

    D.llgrad = @llgrad;
    function [dll, store] = llgrad(theta, data, store)
        
        data = mxe_readdata(data, false);
        index = data.index;
        weight = data.weight;
        datamat = data.data;
        
        if nargin < 3
            store = struct;
        end
        
        store = weighting_intermediate_params(theta, data, store);
 
        data.data = [data.data(1:datadim, index); store.hX];
        data.index =[];
        [dll.G, store] = GateD.llgrad(theta.G, data, store);
        
        % Given the weighting calculate the gradient
        dll.D = cell(num, 1);
        for k = 1:num
            if ~isempty(weight)
                component_weights = store.pcond(k,:) .* weight;
            else
                component_weights = store.pcond(k,:);
            end
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            [dll.D{k}, store.componentStores{k}] = ...
                Experts{k}.llgrad(theta.D{k}, ...
                struct('data',datamat, 'weight',component_weights, 'index',index), ...
                store.componentStores{k});
        end
        
    end

%% |llgraddata|
% See <doc_distribution_common.html#8 distribution structure common members>.

    D.llgraddata = @llgraddata;
    function [dld, store] = llgraddata(theta, data, store)
        
        error('llgraddata not implemented yet');
    end

%% |predict|
%

    D.predict = @predict;
    function [label, store] = predict(theta, data, store)
        data = mxe_readdata(data, false);
        if nargin < 3
            store = struct;
        end
        if isfield(Experts{1},'num')
            datam = data.data(1:datadim,:);
            numE= Experts{1}.num();
            ps = zeros(numE, data.size);
            for k = 1:numE
                datap = [datam; ones(1, data.size)*k];
                [pspart, store] = D.llvec(theta, datap, store);
                store = rmfield(store, 'logpcond');
                ps(k,:) = pspart;
            end
            [notused, label] = max(ps, [], 1);
        end
    end

%% |pdf|
% See <doc_distribution_common.html#10 distribution structure common members>.

%% |sample|
% See <doc_distribution_common.html#11 distribution structure common members>.

    D.sample = @sample;
    function data = sample(theta, n)
        
        error('It does not make sense to sample from gate');
    end

%% |randparam|
% See <doc_distribution_common.html#12 distribution structure common members>.

%% |init|
% See <doc_distribution_common.html#13 distribution structure common members>.

    D.init = @init;
    function theta = init(data, varargin)
        data = mxe_readdata(data);
        label = ceil(rand(1, data.size) * num);
        for k = 1:num
            theta.D{k} = Experts{k}.init(data.data(:,label==k));
        end
        theta.G = GateD.init(data);
    end

%%

%% |estimatedefault|
    D.estimatedefault = @estimatedefault;
    function [varargout] = estimatedefault(varargin)
        [varargout{1:nargout}] = moe_estimatedefault(D, varargin{:}); 
    end

%% |penalizerparam|
% See <doc_distribution_common.html#15 distribution structure common members>.
%

    D.penalizerparam = @penalizerparam;
    function penalizer_theta = penalizerparam(data)
        for k = 1:num
            penalizer_theta.D{k} = Experts{k}.penalizerparam(data);
        end
        penalizer_theta.G = GateD.penalizerparam(data);
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
        
        [costP, store] =  GateD.penalizercost(theta.G, penalizer_theta.G, store);
        for k = 1:num
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            [cost, store.componentStores{k}] = ...
                Experts{k}.penalizercost(theta.D{k}, penalizer_theta.D{k}, ...
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
        
        [gradP.G, store] = GateD.penalizergrad(theta.G, penalizer_theta.G, store);
        for k = 1:num
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            [gradP.D{k}, store.componentStores{k}] = ...
                Experts{k}.penalizergrad(theta.D{k}, penalizer_theta.D{k}, ...
                store.componentStores{k});
        end
       
    end

%% |sumparam|
% See <doc_distribution_common.html#18 distribution structure common members>.

    D.sumparam = @sumparam;
    function theta = sumparam(theta1, theta2)
        theta.D = cell(num, 1);
        for k = 1:num
            theta.D{k} = Experts{k}.sumparam(theta1.D{k}, theta2.D{k});
        end
        theta.G = GateD.sumparam(theta1.G, theta2.G);
    end

%% |scaleparam|
% See <doc_distribution_common.html#19 distribution structure common members>.

    D.scaleparam = @scaleparam;
    function theta = scaleparam(scalar, theta)

        for k = 1:num
            theta.D{k} = Experts{k}.scaleparam(scalar, theta.D{k});
        end
        theta.G = GateD.scaleparam(scalar, theta.G);
    end

%% |sumgrad|
% See <doc_distribution_common.html#20 distribution structure common members>.

%% |scalegrad|
% See <doc_distribution_common.html#21 distribution structure common members>.

%% |entropy|
% See <doc_distribution_common.html#22 distribution structure common members>.

    D.entropy = @entropy;
    function h = entropy(theta)
        error('entropy not implemented yet.');       
    end

%% |kl|
% See <doc_distribution_common.html#23 distribution structure common members>.

    D.kl = @kl;
    function kl = kl(theta)
        error('KL-divergence not implemented yet.');       
    end
    
%% |AICc|
% See <doc_distribution_common.html#24 distribution structure common members>.

%% |BIC|
% See <doc_distribution_common.html#25 distribution structure common members>.

%% |display|
% See <doc_distribution_common.html#26 distribution structure common members>.

    D.display = @display;
    function str = display(theta)
        str = sprintf('Weights (%d-by-%d): %s\n', size(theta.W,1), size(theta.W,2), mat2str(theta.W, 4));
        
        if nargout == 0
            str = [sprintf('%s distribution parameters:\n', D.name()), str];
        end
    end

%%

    D = mxe_addsharedfields(D);
end
