%% |softmaxfactory|
% Construct a softmax gate distribution structure
%
% *Syntax*
%
%   D = softmaxfactory(num, datadim)
%   D = softmaxfactory(num)
%
% *Description*
%
% |D = softmaxfactory(num, datadim)| returns a structure representing a
% softmax gate distribution. |num| is the number of factors. |datadim| is
% the dimensionality of the input space.
%
% |D = softmaxfactory(num)| is the same as above with |datadim = 1|.
%
% *Distribution Parameters*
%
% * *|W|* (|num-by-datadim| matrix) : 
% A matrix containing weights corresponding to each distribution.
%
% *Discrete Probability Density Function Conditioned on Input Space*
%
% The distribution has the following density:
% 
% $$ f(y | x) =  \sum_{k=1}^{num}{ \pi_k(x) \exp(y_k)} $$
%
% where 
%
% $$ \pi_k(y) =  \exp(w_k x) / \sum_{k=1}^{num}{\exp(w_k x)} $$
%
% wherein $num$ is the number of components, $w_i,\ i>1$ are parameters,
% $w_1 = 0$ and $y_k$ is the log-likelihood of the kth conditional
% distribution.
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

function  D = softmaxfactory(num, datadim)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'softmax';

%%

    if nargin < 2
        datadim = 1;
    end

%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.

    D.M = productmanifold(struct('W', euclideanfactory(num-1, datadim+1)));

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
        dim = (num-1) * (datadim+1);
    end

%% |datadim|
% See <doc_distribution_common.html#4 distribution structure common members>.

    D.datadim = @() datadim(); % data space dimensions

%%

    function store = weighting_intermediate_params(theta, data, store)
    % calculate intermediate parameters for ll #1
        
        if ~isfield(store, 'dataTW')
            data = mxe_readdata(data);
            n = data.size;
            data = data.data;
            store.dataTW = theta.W * [data(1:datadim,:); zeros(1,n)];
        end
    end        


    function store = weighting_intermediate_params2(theta, data, store)
        % calculate intermediate parameters for ll #2
        
        
        store = weighting_intermediate_params(theta, data, store);
        
        if ~isfield(store, 'logpgate')
            logpgate = [zeros(1, size(store.dataTW,2)); store.dataTW];
            store.logpgatesum = logsumexp(logpgate);
        end
            
        if ~isfield(store, 'logpcond')
            data = mxe_readdata(data);
            data = data.data;
            store.logpcond = logpgate + data(datadim+1:end,:);
            store.logpcondsum = logsumexp(store.logpcond);
        end
    end

%% |llvec|
% See <doc_distribution_common.html#6 distribution structure common members>.

    D.llvec = @llvec;
    function [llvec, store] = llvec(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        weight = mxe_readweight(data);
        
        store = weighting_intermediate_params2(theta, data, store);
        
        llvec = store.logpcondsum - store.logpgatesum;
        
        if ~isempty(weight)
            llvec = llvec .* weight;
        end
        
    end

%% |ll|
% The same as 
% <doc_distribution_common.html#5 distribution structure common members>
% with the difference that here the log-likelihood of all discrete output
% values is computed as a function of the input vector.
%

    D.ll = @ll;
    function [ll, store] = ll(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        
        [llvec, store] = D.llvec(theta, data, store);
         
        ll = sum(llvec, 2);

    end

%% |llgrad|
% See <doc_distribution_common.html#7 distribution structure common members>.
% Here the gradient of the parameters are given by  
% $$ \sum_{k=1}^{num} a_i log (\frac{f_i( w_i x)}{ 1 + \sum_{k=1}^{num-1} f_k( w_k x)}) $$
%

    D.llgrad = @llgrad;
    function [dll, store] = llgrad(theta, data, store)
        
        data = mxe_readdata(data);
        weight = data.weight;
        n = data.size;
        data = data.data;
        
        if nargin < 3
            store = struct;
        end
        
        store = weighting_intermediate_params2(theta, data, store);
        
        pgate = exp(bsxfun(@minus, store.dataTW, store.logpgatesum));
        
        store.pcond = exp(bsxfun(@minus, store.logpcond, store.logpcondsum));     
        
        part = store.pcond(2:end,:) - pgate;
        
        datap = [data(1:datadim,:); zeros(1,n)];
        
        if ~isempty(weight)
            datap = bsxfun(@times, weight, datap);
        end

        dll.W = part * datap.';
        
    end

%% |llgraddata|
% See <doc_distribution_common.html#8 distribution structure common members>.

    D.llgraddata = @llgraddata;
    function [dld, store] = llgraddata(theta, data, store) %#ok<STOUT,INUSD>
        
        error('llgraddata not implemented yet');
    end
    
% %% |pdf|
% % See <doc_distribution_common.html#10 distribution structure common members>.
% 
% %% |sample|
% % See <doc_distribution_common.html#11 distribution structure common members>.
% 
%     D.sample = @sample;
%     function data = sample(theta, n)
%         
%         error('It does not make sense to sample from gate');
%     end
% 
%% |randparam|
% See <doc_distribution_common.html#12 distribution structure common members>.

%% |init|
% See <doc_distribution_common.html#13 distribution structure common members>.

    D.init = @init;
    function theta = init(data, varargin) %#ok<INUSD>
        data = mxe_readdata(data);
        
        optionsinit.verbosity = 2;
        optionsinit.maxiter = 50;
        optionsinit.solver = 'cg';
        optionsinit.penalize = true;
        optionsinit.tolcostdiff = 1e-5;
        
        label = ceil(rand(1, data.size) * num);
        
        Dsmax = mnlfactory(datadim, num);
        optionsinit.regularize = false;
        theta= Dsmax.estimate([data.data(1:datadim,:); label], optionsinit);
        %theta = D.M.rand();
    end

%%

% %% |estimatedefault|
%     D.estimatedefault = @estimatedefault;
%     function [varargout] = estimatedefault(varargin)
%         [varargout{1:nargout}] = factorial_estimatedefault(D, varargin{:}); 
%     end

%% |estimateMstep|
% When using EM algorithm for MoE, M-step estimation for gates is needed
% The cost function here is convex and therefore it converges rapidly
    D.estimateMstep = @estimateMstep;
    function [varargout] = estimateMstep(varargin)
        [varargout{1:nargout}] = softmax_estimateMstep(D, varargin{:}); 
    end

%% |penalizerparam|
% See <doc_distribution_common.html#15 distribution structure common members>.

    D.penalizerparam = @penalizerparam;
    function penalizer_theta = penalizerparam(data) %#ok<INUSD>
        penalizer_theta = [];
    end

%% |penalizercost|
% See <doc_distribution_common.html#16 distribution structure common members>.

    D.penalizercost = @penalizercost;
    function [costP, store] = penalizercost(theta, penalizer_theta, store) %#ok<INUSL>
        
        costP = 0;
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
        
        gradP.W = zeros(size(theta.W));
       
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
