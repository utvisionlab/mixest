%% |mvnicholfactory|
% Construct a multi-variate normal distribution structure
%
% *Syntax*
%
%   D = mvnicholfactory(datadim)
%
% *Description*
%
% |D = mvnicholfactory(datadim)| returns a structure representing a
% |datadim|-dimensional normal distribution.
%
% *Distribution Parameters*
%
% * *|mu|* (|datadim-by-1| vector) : The mean vector.
% * *|cholsinv|* (|datadim-by-datadim| matrix) : Cholesky factor of inverse
% covariance matrix (or the variance, $\sigma^2$, when |datadim|=1).
%
% *Probability Density Function*
%
% The distribution has the following density:
% 
% $$ f(x;\mu,\Sigma)=
% (2\pi)^{-\frac{n}{2}}|\Sigma|^{-\frac{1}{2}}
% \exp\left(-\frac{1}{2}({x}-{\mu})^T{\Sigma}^{-1}({x}-{\mu})
% \right) $$
%
% where $n$ is the data space dimensions, $\mu$ is the mean vector and
% $\Sigma$ is the covariance matrix.
% 
% *Example*
%
%   % Construct a bivariate normal distribution:
%   D = mvnicholfactory(2);
%   % Build a parameter structure for it:
%   theta = struct('mu', [0; 0], 'sigma', [sqrt(0.5) 0; 0 sqrt(0.5)]);
%   % Plot the PDF:
%   x = -5:0.2:5;
%   y = -5:0.2:5;
%   [X, Y] = meshgrid(x, y);
%   data = [X(:) Y(:)]';
%   f = D.pdf(theta, data);
%   surf(X, Y, reshape(f, size(X)));
%
% <<img/mvnfactory_01.png>>
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

function D = mvnicholfactory(datadim)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'mvnichol';

%%

    assert(datadim >= 1, 'datadim must be an integer larger than or equal to 1.');
    Dmvn = mvnfactory(datadim);

%% 
% Flag to control the memory usage (resulting code will be slower)

    low_memory = true;

%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.

    muM = euclideanfactory(datadim);
    cholsinvM = triufactory(datadim, datadim);
    D.M = productmanifold(struct('mu', muM, 'cholsinv', cholsinvM));

%% |dim|
% See <doc_distribution_common.html#3 distribution structure common members>.

    D.dim = @() datadim*(datadim+1)/2 + datadim; % parameter space dimensions

%% |datadim|
% See <doc_distribution_common.html#4 distribution structure common members>.

    D.datadim = @() datadim; % data space dimensions

%% |main2new|
% Convert usual MVN parameters to reparameterized parameters
%
% *Syntax*
%
%   theta = D.main2new(theta)
%

    D.main2new = @main2new;
    function theta = main2new(theta)
        if isfield(theta,'mu')
            % computing the cholesky factor of the inverse covariance
            R = chol(theta.sigma); % C=R' R ( R upper trangular)
            theta.cholsinv = R \ eye(datadim);
            theta = rmfield(theta,'sigma');
        end
    end

%% |new2main|
% Convert reparameterized parameters to usual MVN parameters
%
% *Syntax*
%
%   theta = D.new2main(theta)
%

    D.new2main = @new2main;
    function theta = new2main(theta)
        if isfield(theta,'cholsinv');
            R = theta.cholsinv \ eye(datadim);
            theta.sigma = R.' * R;
            theta = rmfield(theta,'cholsinv');
        end
    end

%%

    function store = ll_intermediate_params2(theta, data, store)
    % calculate intermediate parameters for ll #2
    % Note: store must contain fields created by the previous ll_intermediate_params functions

        if ~isfield(store, 'data') %|| ~isfield(store, 'weight')
            data = mxe_readdata(data);
            data = data.data;
            
            % Subtracting mean from the data
            data = bsxfun(@plus, data , -theta.mu);
            store.data = data;
        end
    end
    function store = ll_intermediate_params3(weight, store)
    % calculate intermediate parameters for ll #3
    % Note: store must contain fields created by the previous ll_intermediate_params functions
        
        if ~isfield(store, 'sumW') || ~isfield(store, 'DDT')
            data = store.data;
            
            % Multiply data by the square root of the weight
            if ~isempty(weight)
                store.sumW = sum(weight, 2);
                weight = sqrt(weight);
                % Rinv'*data is really slow, so it is better to do the following
                data = bsxfun(@times, data, weight);
            else
                store.sumW = size(data, 2);
            end
            
            % Calculate weight*data*data'
            store.DDT = (data * data.');
            
            store.data = data;
            store.weight = weight;
        end
    end

%% |ll|
% See <doc_distribution_common.html#5 distribution structure common members>.

    D.ll = @ll;
    function [ll, store] = ll(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        weight = mxe_readweight(data);

        store = ll_intermediate_params2(theta, data, store);
        logdetR = sum(log(abs(diag(theta.cholsinv))));
        if ~isfield(store,'Sinv')
            store.Sinv = theta.cholsinv * theta.cholsinv.';
        end
        Sinv = store.Sinv;

        store = ll_intermediate_params3(weight, store);
        sumW = store.sumW;
        DDT = store.DDT;
        
        % Calculate u = sum_i(w_i d_i' C^-1 d_i)
        u = DDT(:).'*Sinv(:);

        % Calculate the log-likelihood
        ll = sumW *(- (0.5*datadim)*log(2*pi) + logdetR) - (1/2)*u;
        
    end

%% |llvec|
% See <doc_distribution_common.html#6 distribution structure common members>.

    D.llvec = @llvec;
    function [llvec, store] = llvec(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        weight = mxe_readweight(data);

        store = ll_intermediate_params2(theta, data, store);
        logdetR = sum(log(abs(diag(theta.cholsinv))));
        
        % u = X' C^-1 X 
        u = sum((theta.cholsinv' * store.data).^2, 1);
        
        if low_memory % if size of data is large, better to remove data in store
            store = rmfield(store,'data');
        end
        
        llvec = - (0.5*datadim)*log(2*pi) + logdetR - (1/2)*u;            
        if ~isempty(weight)
            llvec = weight .* llvec;
        end
    end

%% |llgrad|
% See <doc_distribution_common.html#7 distribution structure common members>.

    D.llgrad = @llgrad;
    function [dll, store] = llgrad(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        data = mxe_readdata(data);
        weight = data.weight;
        data = data.data;

        store = ll_intermediate_params2(theta, data, store);
        if ~isfield(store,'Sinv')
            store.Sinv = theta.cholsinv * theta.cholsinv.';
        end
        Sinv = store.Sinv;

        store = ll_intermediate_params3(weight, store);
        sumW = store.sumW;
        DDT = store.DDT;

        % gradient with respect to mu
        if isfield(store, 'dld')
            dll.mu = - sum(store.dld, 2);
        else
            if ~isempty(weight)
                data = bsxfun(@times, store.data, store.weight);
            else
                data = store.data;
            end
            dll.mu = Sinv * sum(data, 2);
        end

        % gradient with respect to sigma
        dll.cholsinv = diag(sumW ./diag(theta.cholsinv)) - theta.cholsinv.' * DDT ;
        dll.cholsinv = tril(dll.cholsinv).';
        
        if low_memory % if size of data is large, better to remove data from store
            store = rmfield(store,'data');
            store = rmfield(store,'weight');
        end
    end

%% |llgraddata|
% See <doc_distribution_common.html#8 distribution structure common members>.

    D.llgraddata = @llgraddata;
    function [dld, store] = llgraddata(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        weight = mxe_readweight(data);

        store = ll_intermediate_params2(theta, data, store);
        if ~isfield(store,'Sinv')
            store.Sinv = theta.cholsinv * theta.cholsinv.';
        end
        Sinv = store.Sinv;

        data = store.data;

        if ~isempty(weight)
            data = bsxfun(@times, data, weight);
        end

        % gradient with respect to data
        dld = - Sinv * data;
        
        store.dld = dld;
        
        if low_memory % if size of data is large, better to remove data from store
            store = rmfield(store,'data');
        end
    end
        
%% |cdf|
% See <doc_distribution_common.html#9 distribution structure common members>.

    D.cdf = @cdf;
    function y = cdf(theta, data)
        data = mxe_readdata(data);
        data = data.data;
        if datadim == 1
            y = 0.5 * erfc( (theta.mu - data) ./ sqrt(2*theta.sigma) );
        else
            y = mvncdf(data.', theta.mu.', theta.sigma).'; %TODO fro Genz
        end
    end

%% |pdf|
% See <doc_distribution_common.html#10 distribution structure common members>.

%% |sample|
% See <doc_distribution_common.html#11 distribution structure common members>.

    D.sample = @sample;
    function data = sample(theta, n)
        
        if nargin < 2, n = 1; end
        
        data = randn(datadim, n);
        L = theta.cholsinv / eye(datadim);
        data = L * data; 
        data = bsxfun(@plus, data , theta.mu);
    end

%% |randparam|
% See <doc_distribution_common.html#12 distribution structure common members>.

%% |init|
% See <doc_distribution_common.html#13 distribution structure common members>.

    D.init = @init;
    function theta = init(data, varargin)
        
        theta = estimatedefault(data);
    end

%% |estimatedefault|
% Default estimation function for multi-variate normal distribution. This
% function implements the maximum likelihood method.
%
% *Syntax*
%
%   theta = D.estimatedefault(data)
%   theta = D.estimatedefault(data, options)
%   [theta, D] = D.estimatedefault(...)
%   [theta, D, info] = D.estimatedefault(...)
%   [theta, D, info, options] = D.estimatedefault(...)
%
% *Description*
%
% |theta = D.estimatedefault(data)| returns estimated parameters for the
% distribution |D|, using |data|.
%
% |theta = D.estimatedefault(data, options)| utilizes applicable options
% from the |options| structure in the estimation procedure.
%
% |[theta, D] = D.estimatedefault(...)| also returns |D|, the distribution
% structure for which |theta| is applicable. (This is the same as the
% distribution structure |D| from which you called |estimate|, and so it
% should not normally be used. The purpose of including it in the output is
% to maintain compatibility with other estimation functions).
%
% |[theta, D, info] = D.estimatedefault(...)| also returns |info|, a
% structure array containing information about successive iterations
% performed by iterative estimation functions.
%
% |[theta, D, info, options] = D.estimatedefault(...)| also returns the
% effective |options| used, so you can see what default values the function
% used on top of the options you possibly specified.
%
% For information about the output |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>. You may also want to read about
% <../estimation_options.html |options|> or
% <../estimation_statistics_structure.html |info|> arguments.
%
% *Available Options*
%
% Currently no options are available for this function.
%
% *Returned |info| fields*
%
% The method used is not iterative and so the returned |info| is empty.
%
% *Example*
%
%   % create a Gaussian distribution
%   D = mvnfactory(1);
%   % generate 1000 random data points
%   data = randn(1,1000) .* 2 + 1;
%   % estimate distribution parameters to fit the data
%   theta = D.estimatedefault(data)
%

    D.estimatedefault = @estimatedefault;
    function [varargout] = estimatedefault(varargin)
        [varargout{1:nargout}] = mvn_estimatedefault(Dmvn, varargin{:});
        varargout{1} = main2new(varargout{1});
    end

%% |penalizerparam|
% See <doc_distribution_common.html#15 distribution structure common members>.



%% |penalizercost|
% See <doc_distribution_common.html#16 distribution structure common members>.


%% |penalizergrad|
% See <doc_distribution_common.html#17 distribution structure common members>.


%% |sumparam|
% See <doc_distribution_common.html#18 distribution structure common members>.

    D.sumparam = @sumparam;
    function theta = sumparam(theta1, theta2)
        theta.mu = theta1.mu + theta2.mu;
        theta.cholsinv = theta1.cholsinv + theta2.cholsinv;
    end

%% |scaleparam|
% See <doc_distribution_common.html#19 distribution structure common members>.

    D.scaleparam = @scaleparam;
    function theta = scaleparam(scalar, theta)
        theta.mu = scalar * theta.mu;
        theta.cholsinv = scalar * theta.cholsinv;
    end

%% |sumgrad|
% See <doc_distribution_common.html#20 distribution structure common members>.

%% |scalegrad|
% See <doc_distribution_common.html#21 distribution structure common members>.

%% |entropy|
% See <doc_distribution_common.html#22 distribution structure common members>.

    D.entropy = @entropy;
    function h = entropy(theta)
        logdetR = -sum(log( diag( theta.cholsinv ) )); 
        h = logdetR  + 0.5 * datadim * ( log(2*pi)  + 1 );            
    end

%% |kl|
% See <doc_distribution_common.html#23 distribution structure common members>.

    D.kl = @kl;
    function kl = kl(theta1, theta2)
        % KL-divergance is given by 
        % 1/2*[ trace(M) + (mu1-mu2)^T M (mu1-mu2) - n - log(det(M))  ]
        % where M = Sigma_2^-1 *Sigma1 (index 1 is true and 2 is model)
        
        % A fast computation of both M and log(det(M)) using cholesky
        U = theta1.cholsinv; % Y.sigma = U'*U (U is upper-traingular)
        V = theta2.cholsinv;
        U = U \ eye(datadim); % V = inv_triu(V);
        logdetR = sum(log( diag( U ) )) + sum(log( diag( V ) ));
        mu_diff = theta2.mu - theta1.mu; 
        term = sum((mu_diff'*V).^2); %mu_diff'*inv(sigma_2)*mu_diff
        V =  U * V; % trace(AA')=sum(A(:).^2)
        trM = sum(V(:).^2); % trM = trace(Sigma_2^-1 *Sigma1)
        kl = 0.5*(trM + term - datadim) - logdetR;
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
        str = [str, sprintf('mean (%d-by-1): %s\n', size(theta.mu,1), mat2str(theta.mu, 4))];
        str = [str, sprintf('cholsinv (%d-by-%d): %s\n', size(theta.cholsinv,1), size(theta.cholsinv,2), mat2str(theta.cholsinv, 4))];
        
        if nargout == 0
            str = [sprintf('%s distribution parameters:\n', D.name()), str];
        end
    end

%% |selfsplit|
% See <doc_distribution_common.html#27 distribution structure common members>.

    
%% |selfmerge|
% See <doc_distribution_common.html#28 distribution structure common members>.


%% |visualize|
%
% *Syntax*
%
%   handle_array = D.visualize(D, theta, vis_options)
%

    D.visualize = @visualize;
    function [varargout] = visualize(varargin)
        varargin{1} = new2main(varargin{1});
        [varargout{1:nargout}] = mvn_visualize(D, varargin{:});
    end

%%

    D = mxe_addsharedfields(D);
end
