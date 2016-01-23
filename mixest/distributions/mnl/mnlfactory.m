%% |mnlfactory|
% Construct a multinomial logistic distribution (MaxEnt classifier) 
%
% *Syntax*
%
%   D = mnlfactory(datadim, num)
%   D = mnlfactory(datadim)
%
% *Description*
%
% |D = mnlfactory(datadim, num)| returns a structure representing a
% multi-nomial logit distribution. |datadim| is the dimension of input
% space and |num| is the number of labels.
%
% |D = mnlfactory(datadim)| is the same as above with |num = 2|, this case
% equals usual logistic regression.
%
% *Distribution Parameters*
%
% * *|W|* (|(num - 1)-by-datadim| matrix) : 
% A matrix containing weights from input space to |num| classifier
%
% *Probability Density Function*
%
% The distribution has the following density:
% 
% $$ f(l | x)= w_{l-1}^T x / (1 + \sum_{k=2}^{num}{ \exp(w_{k-1}^T x)}) $$
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

function D = mnlfactory(datadim, num)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'multinomial logit';
    
    if nargin < 2
        num = 2;
    end
    
%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.

    D.M = productmanifold(struct('W', euclideanfactory(num-1, datadim+1)));

%% |num|
% Number of components
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
        dim = (num-1) * (datadim+1);
    end

%% |datadim|
% See <doc_distribution_common.html#4 distribution structure common members>.

    D.datadim = @() datadim(); % data space dimensions


%%     
    function store = weighting_intermediate_params(theta, data, store)
        % calculate intermediate parameters for ll #2
        
        if ~isfield(store, 'logpsum')
            data = mxe_readdata(data);
            
            store.dataTW = theta.W * [data.data(1:datadim,:);ones(1,data.size)];
            logp = [zeros(1, data.size); store.dataTW];
            store.logpsum = logsumexp(logp);
            
        end
    end

%% |llvec|
% See <doc_distribution_common.html#6 distribution structure common members>.

    D.llvec = @llvec;
    function [llvec, store] = llvec(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        data = mxe_readdata(data);
        
        store = weighting_intermediate_params(theta, data, store);
        label = data.data(end, :);
        labg1ind = label > 1;
        label(labg1ind) = label(labg1ind) - 1;
        index = sub2ind([num-1 data.size], label, 1:data.size);
        llvec = zeros(1, data.size);
        llvec(labg1ind) = store.dataTW(index(labg1ind)) - ...
            store.logpsum(labg1ind);
        llvec(~labg1ind) = - store.logpsum(~labg1ind);
        
        if ~isempty(data.weight)
            llvec = llvec .* data.weight;
        end
        
    end

%% |ll|
% See <doc_distribution_common.html#5 distribution structure common members>.

    D.ll = @ll;
    function [ll, store] = ll(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        data = mxe_readdata(data);
        
        if isempty(data.weight)
            store = weighting_intermediate_params(theta, data, store);
            label = data.data(end, :);
            labg1ind = label > 1;
            label(labg1ind) = label(labg1ind) - 1;
            index = sub2ind([num-1 data.size], label, 1:data.size);
            ll =  sum(store.dataTW(index(labg1ind))) - sum(store.logpsum);
        else
            [llvec, store] = D.llvec(theta, data, store);
            ll = sum(llvec, 2);
        end
    end

%% |llgrad|
% See <doc_distribution_common.html#7 distribution structure common members>.

    D.llgrad = @llgrad;
    function [dll, store] = llgrad(theta, data, store)
        
        data = mxe_readdata(data);
        
        if nargin < 3
            store = struct;
        end
        
        store = weighting_intermediate_params(theta, data.data, store);
        
        label = data.data(end, :);
        datap = [data.data(1:end-1, :); ones(1,data.size)];
        if ~isempty(data.weight)
            datap = bsxfun(@times, data.weight, datap);
        end
        part = exp(bsxfun(@minus, store.dataTW, store.logpsum));
        
        dll.W = part * datap.';
        
        W = zeros(num-1, datadim+1);
        for k = 2:num
            index = (label == k);
            W(k-1,:) = sum(datap(:,index), 2).';
        end
        dll.W = W - dll.W;
        
    end

%% |llgraddata|
% See <doc_distribution_common.html#8 distribution structure common members>.

    D.llgraddata = @llgraddata;
    function [dld, store] = llgraddata(theta, data, store) %#ok<STOUT,INUSD>
        
        error('llgraddata not implemented yet');
    end

%% |predict|
%

    D.predict = @predict;
    function [label, store] = predict(theta, data, store)
        data = mxe_readdata(data, false);
        datam = data.data(1:datadim,:);
        ps = zeros(num, data.size);
        for k = 1:num
            datap = [datam; ones(1, data.size)*k];
            [pspart, store] = D.llvec(theta, datap);
            ps(k,:) = pspart;
        end
        [notused, label] = max(ps, [], 1);
    end

%% |pdf|
% See <doc_distribution_common.html#10 distribution structure common members>.

%% |sample|
% See <doc_distribution_common.html#11 distribution structure common members>.

    D.sample = @sample;
    function data = sample(theta, n) %#ok<STOUT,INUSD>
        
        error('sampling not implemented yet');
    end

%% |randparam|
% See <doc_distribution_common.html#12 distribution structure common members>.

%% |init|
% See <doc_distribution_common.html#13 distribution structure common members>.

    D.init = @init;
    function theta = init(data, varargin) %#ok<INUSD>
        theta = D.M.rand();
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

    D.penalizerparam = @penalizerparam;
    function penalizer_theta = penalizerparam(data) %#ok<INUSD>
        penalizer_theta = 1;
    end

%% |penalizercost|
% See <doc_distribution_common.html#16 distribution structure common members>.

    D.penalizercost = @penalizercost;
    function [costP, store] = penalizercost(theta, penalizer_theta, store) %#ok<INUSL>
        
        costP = - 0.5 * penalizer_theta * norm(theta.W(:,1:end-1), 'fro')^2;
        if nargin < 3
            store = struct;
        end
     
    end

%% |penalizergrad|
% See <doc_distribution_common.html#17 distribution structure common members>.

    D.penalizergrad = @penalizergrad;
    function [gradP, store] = penalizergrad(theta, penalizer_theta, store) %#ok<INUSL>
        
        if nargin < 3
            store = struct;
        end
        gradP.W = - penalizer_theta * theta.W;
        gradP.W(:,end) = 0;
       
    end

%% |sumparam|
% See <doc_distribution_common.html#18 distribution structure common members>.

    D.sumparam = @sumparam;
    function theta = sumparam(theta1, theta2)

        theta.W = theta1.W + theta2.W;
    end

%% |scaleparam|
% See <doc_distribution_common.html#19 distribution structure common members>.

    D.scaleparam = @scaleparam;
    function theta = scaleparam(scalar, theta)

        theta.W = scalar * theta.W;
    end

%% |sumgrad|
% See <doc_distribution_common.html#20 distribution structure common members>.

%% |scalegrad|
% See <doc_distribution_common.html#21 distribution structure common members>.

%% |entropy|
% See <doc_distribution_common.html#22 distribution structure common members>.

    D.entropy = @entropy;
    function h = entropy(theta) %#ok<STOUT,INUSD>
        error('entropy not implemented yet.');       
    end

%% |kl|
% See <doc_distribution_common.html#23 distribution structure common members>.

    D.kl = @kl;
    function kl = kl(theta) %#ok<STOUT,INUSD>
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