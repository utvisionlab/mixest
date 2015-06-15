%% |cmvnfactory|
% Construct a conditional multi-variate normal distribution structure
%
% *Syntax*
%
%   D = cmvnfactory(datadim, conddim)
%
% *Description*
%
% |D = cmvnfactory(datadim, conddim)| returns a structure representing a
% conditional normal distribution. |datadim| is the data dimensions and
% |conddim| is the condition dimensions.
%
% *Distribution Parameters*
%
% * *|W|* (|(conddim)-by-(datadim+1)| matrix) : 
% A matrix containing weights from input space to output.
% * *|sigma|* (|conddim-by-conddim| matrix) : The covariance matrix (or the
% variance, $\sigma^2$, when |datadim|=1).
%
% *Probability Density Function*
%
% The distribution has the following density:
% 
% $$ f(y|x;\mu,\Sigma)=
% (2\pi)^{-\frac{n}{2}}|\Sigma|^{-\frac{1}{2}}
% \exp\left(-\frac{1}{2}({y}-{\mu_x})^T{\Sigma}^{-1}({y}-{\mu_x})
% \right) $$
%
% where $n$ is the data space dimensions, $\Sigma$ is the covariance
% matrix and $\mu_x$ is the mean vector computed by
%
% $$ \mu_x = W [x;1] $$
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

function D = cmvnfactory(datadim, conddim)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'cmvn';

%%

    assert(datadim >= 1, 'datadim must be an integer larger than or equal to 1.');

%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.

    W = euclideanfactory(conddim, datadim+1);
    sigmaM = spdfactory(conddim);
    D.M = productmanifold(struct('W', W, 'sigma', sigmaM));

%% |dim|
% See <doc_distribution_common.html#3 distribution structure common members>.

    D.dim = @() conddim*(conddim+1)/2 + conddim*(datadim +1); % parameter space dimensions

%% |datadim|
% See <doc_distribution_common.html#4 distribution structure common members>.

    D.datadim = @() datadim; % marginal data space dimensions
    
%% |conddim|
%

    D.conddim = @() conddim; % data space dimension

%%

    function store = ll_intermediate_params1(theta, store)
    % calculate intermediate parameters for ll #1 (parameters depending only on theta)
        
        if ~isfield(store, 'Rinv') || ~isfield(store, 'logdetR')
            % Compute inverse cholesky and 1/2 log(det(.)) of it
            [R, p] = chol(theta.sigma); % C=R' R ( R upper trangular)
            if p > 0
                warning('matrix is not symmetric positive definite ...');
                t = 1e-10;
                while p > 0
                    theta.sigma = theta.sigma + t * eye(conddim);
                    [R, p] = chol(theta.sigma); %#ok<ASGLU>
                    t = t * 100;
                end
                R = chol(theta.sigma);
            end
            store.Rinv = R \ eye(conddim); % Faster version of Rinv = inv_triu(R);
            store.logdetR = sum(log(diag(R)));
        end
        if ~isfield(store, 'Sinv')
            Rinv = store.Rinv;
            % Non-vectorized output
            store.Sinv = Rinv * Rinv';
        end
    end
    function store = ll_intermediate_params2(theta, data, store)
    % calculate intermediate parameters for ll #2
    % Note: store must contain fields created by the previous ll_intermediate_params functions

        if ~isfield(store, 'data') %|| ~isfield(store, 'weight')
            data = mxe_readdata(data);
            n = data.size;
            data = data.data;
            
            % Subtracting mean from the data
            store.mu = theta.W * [data(1:datadim,:); ones(1, n)];
            data = bsxfun(@plus, data(datadim+1:end,:) , -store.mu);
            store.data = data;
            %store.weight = weight;
        end
    end
    function store = ll_intermediate_params3(weight, store)
    % calculate intermediate parameters for ll #3
    % Note: store must contain fields created by the previous ll_intermediate_params functions
        
        if ~isfield(store, 'sumW') || ~isfield(store, 'DDT')
            data = store.data;
            %weight = store.weight;
            
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
        
        store = ll_intermediate_params1(theta, store);
        logdetR = store.logdetR;
        Sinv = store.Sinv;

        store = ll_intermediate_params2(theta, data, store);

        store = ll_intermediate_params3(weight, store);
        sumW = store.sumW;
        DDT = store.DDT;
        
        % Calculate u = sum_i(w_i d_i' C^-1 d_i)
        u = DDT(:).'*Sinv(:);

        % Calculate the log-likelihood
        ll = sumW *(- (0.5*conddim)*log(2*pi) - logdetR) - (1/2)*u;
    end

%% |llvec|
% See <doc_distribution_common.html#6 distribution structure common members>.

    D.llvec = @llvec;
    function [llvec, store] = llvec(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        weight = mxe_readweight(data);
        
        store = ll_intermediate_params1(theta, store);
        Rinv = store.Rinv;
        logdetR = store.logdetR;

        store = ll_intermediate_params2(theta, data, store);
        
        % u = X' C^-1 X 
        u = sum((Rinv' * store.data).^2, 1);            
        llvec = - (0.5*conddim)*log(2*pi) - logdetR - (1/2)*u;            
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
        datamat = [data.data(1:datadim,:); ones(1, data.size)];
        data = data.data;
        
        store = ll_intermediate_params1(theta, store);
        Sinv = store.Sinv;

        store = ll_intermediate_params2(theta, data, store);

        store = ll_intermediate_params3(weight, store);
        sumW = store.sumW;
        DDT = store.DDT;

        % gradient with respect to mu
        if ~isfield(store, 'dld')
            if ~isempty(weight)
                data = bsxfun(@times, store.data, store.weight);
            else
                data = store.data;
            end
            store.dld = -Sinv * data;
        end
        dll.W = -store.dld * datamat.';
        % gradient with respect to sigma
        dll.sigma = 1/2 * ( -sumW * Sinv + Sinv * DDT *Sinv);
    end

%% |llgraddata|
% See <doc_distribution_common.html#8 distribution structure common members>.

    D.llgraddata = @llgraddata;
    function [dld, store] = llgraddata(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        weight = mxe_readweight(data);
        
        store = ll_intermediate_params1(theta, store);
        Sinv = store.Sinv;

        store = ll_intermediate_params2(theta, data, store);

        data = store.data;

        if ~isempty(weight)
            data = bsxfun(@times, data, weight);
        end

        % gradient with respect to data
        dld = - Sinv * data;
        
        store.dld = dld;
    end
        
%% |cdf|
% See <doc_distribution_common.html#9 distribution structure common members>.

    D.cdf = @cdf;
    function y = cdf(theta, data)
        data = mxe_readdata(data);
        data = data.data;
        mu = theta.W * [data(1:datadim,:); zeros(1,n)];
        y = mvncdf(data(datadim+1:end,:).', mu.', theta.sigma).'; %TODO
    end

%% |pdf|
% See <doc_distribution_common.html#10 distribution structure common members>.

%% |sample|
% See <doc_distribution_common.html#11 distribution structure common members>.

    D.sample = @sample;
    function data = sample(theta, datamargin)
        
        if isscalar(datamargin)
            n = datamargin;
            datamargin = randn(datadim, n);
        end
        n = size(datamargin,2);
        data = randn(conddim, n);
        L = chol(theta.sigma,'lower');
        data = L * data;
        mu = theta.W * [datamargin; zeros(1,n)];
        data = [datamargin; bsxfun(@plus, data , mu)];
    end

%% |randparam|
% See <doc_distribution_common.html#12 distribution structure common members>.

%% |init|
% See <doc_distribution_common.html#13 distribution structure common members>.

    D.init = @init;
    function theta = init(data, varargin)
        initoptions = mxe_options();
        theta = estimatedefault(data, initoptions);
    end

%% |estimatedefault|
% Default estimation function for conditional multi-variate normal distribution. This
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
%   % create a cmvn distribution
%   D = cmvnfactory(1, 1);
%   % generate 1000 random data points
%   data = randn(2,1000) .* 2 + 1;
%   % estimate distribution parameters to fit the data
%   theta = D.estimatedefault(data)
%

    D.estimatedefault = @estimatedefault;
    function [theta, Dout, info, options] = estimatedefault(data, options)
        
        if nargin < 2
            options = mxe_options();
        else
            options = mxe_options(options);
        end
        data = mxe_readdata(data);
        weight = data.weight;
        n = data.size;
        if options.penalize
            penalizer_est_theta = penalizerparam(data);
        end
        data = [data.data(1:datadim,:); ones(1, n); data.data(datadim+1:end,:)];
        if ~isempty(weight)
            data = bsxfun(@times, sqrt(weight), data);
        end
        sigma = (data * data.')/n;
        sigma1 = sigma(datadim+2:end,datadim+2:end );
        sigma2 = sigma(1:datadim+1, 1:datadim+1);
        sigma12 = sigma(datadim+2:end, 1:datadim+1);
        theta.W = sigma12 / sigma2;
        theta.sigma = sigma1 - theta.W * sigma12.';
        if options.penalize
            denum = (n + penalizer_est_theta.nu);
            if isscalar(penalizer_est_theta.invLambda)
                 theta.sigma = n/ denum * theta.sigma + ...
                     1/denum * diag( repmat(penalizer_est_theta.invLambda, [1 size(theta.sigma,1)]));
            else
                theta.sigma = n/ denum * theta.sigma + ...
                    1/denum * penalizer_est_theta.invLambda;
            end
        end
        Dout = D;
        if nargout > 2
            info = [];
            if nargout > 3
                options = struct;
            end
        end
    end

%% |penalizerparam|
% See <doc_distribution_common.html#15 distribution structure common members>.
%
% *Penalizer Info*
%
% The default penalizer for this distribution is the Inverse-Wishart
% distribution for covariance and Normal distribution for mean
%
% Inverse-Wishart prior on covariance has the following form:
%
% $$ f(\Sigma) = 
% |\Sigma|^{-(nu+d+1)/2} \exp(-0.5 trace(\Sigma^{-1} invLambda)) $$
%
% where
%
% * *|nu|* (scalar) : Degrees of freedom.
% * *|invLambda|* (|datadim-by-datadim| matrix) : The inverse scale matrix.
%
% Normal prior on the mean has the following form:
%
% $$ f(\mu|\Sigma) = 
% |\Sigma|^{-1/2} exp(-kappa/2 (\mu-\mu_p)^T \Sigma^{-1} (\mu-\mu_p)) $$
%
% where
%
% * *|mu_p|* (|datadim-by-1| vector) : The mean vector.
% * *|kappa|* (scalar) : the shrinkage parameter
%

    D.penalizerparam = @penalizerparam;
    function penalizer_theta = penalizerparam(data)
        % Using Fraley&Raftery (2007) method for computing parameters
        data = mxe_readdata(data);
        n = data.size;
        data = data.data(datadim+1:end, :);

        penalizer_theta.nu = 2;
        mu = sum(data, 2)/n;
        data = bsxfun(@minus, data, mu);
        sigma = (data * data.')/n;
        % Check if it is multiplication of identity
        if isequal(sigma, sigma(1,1) * eye(conddim))
            penalizer_theta.invLambda = 1/sigma(1,1);
        else
            penalizer_theta.invLambda = inv(sigma);
        end
    end

%% |penalizercost|
% See <doc_distribution_common.html#16 distribution structure common members>.

    D.penalizercost = @penalizercost;
    function [costP, store] = penalizercost(theta, penalizer_theta, store)
        
        if nargin < 3
            store = struct;
        end
        
        store = ll_intermediate_params1(theta, store);
        logdetR = store.logdetR;
        Sinv = store.Sinv;
        
        costP = - (penalizer_theta.nu + datadim + 1) * logdetR;
        % If data to penalizeparam is whitened then invLambda is scalar
        if isscalar(penalizer_theta.invLambda)
            costP = costP - 0.5 * penalizer_theta.invLambda * trace(Sinv);
        else
            costP = costP - 0.5 * (penalizer_theta.invLambda(:).' * Sinv(:));
        end
    end

%% |penalizergrad|
% See <doc_distribution_common.html#17 distribution structure common members>.

    D.penalizergrad = @penalizergrad;
    function [gradP, store] = penalizergrad(theta, penalizer_theta, store)
        
        if nargin < 3
            store = struct;
        end
        
        store = ll_intermediate_params1(theta, store);
        Sinv = store.Sinv;
        
        if isscalar(penalizer_theta.invLambda)
            gradP.sigma = -0.5 * (penalizer_theta.nu + datadim + 1) * Sinv +  ...
                (0.5 * penalizer_theta.invLambda) * (Sinv * Sinv.');
        else
            gradP.sigma = -0.5 * (penalizer_theta.nu + datadim + 1) * Sinv +  ...
                0.5 * Sinv * penalizer_theta.invLambda * Sinv;
        end
        gradP.W = zeros(conddim, datadim+1);
    end

%% |sumparam|
% See <doc_distribution_common.html#18 distribution structure common members>.

    D.sumparam = @sumparam;
    function theta = sumparam(theta1, theta2)
        theta.W = theta1.W + theta2.W;
        theta.sigma = theta1.sigma + theta2.sigma;
    end

%% |scaleparam|
% See <doc_distribution_common.html#19 distribution structure common members>.

    D.scaleparam = @scaleparam;
    function theta = scaleparam(scalar, theta)
        theta.W = scalar * theta.W;
        theta.sigma = scalar * theta.sigma;
    end

%% |sumgrad|
% See <doc_distribution_common.html#20 distribution structure common members>.

%% |scalegrad|
% See <doc_distribution_common.html#21 distribution structure common members>.

%% |entropy|
% See <doc_distribution_common.html#22 distribution structure common members>.

    D.entropy = @entropy;
    function h = entropy(theta)
        logdetR = sum(log( diag( chol(theta.sigma) ) )); 
        h = logdetR  + 0.5 * conddim * ( log(2*pi)  + 1 );            
    end

%% |kl|
% See <doc_distribution_common.html#23 distribution structure common members>.

    D.kl = @kl;
    function kl = kl(theta1, theta2)
        % KL-divergance is given by 
        % 1/2*[ trace(M) + (mu1-mu2)^T M (mu1-mu2) - n - log(det(M))  ]
        % where M = Sigma_2^-1 *Sigma1 (index 1 is true and 2 is model)
        
        % A fast computation of both M and log(det(M)) using cholesky
        U = chol(theta1.sigma); % Y.sigma = U'*U (U is upper-traingular)
        V = chol(theta2.sigma);
        V = V \ eye(conddim); % V = inv_triu(V);
        logdetR = sum(log( diag( U ) )) + sum(log( diag( V ) ));
        mu_diff = theta2.mu - theta1.mu; 
        term = sum((mu_diff'*V).^2); %mu_diff'*inv(sigma_2)*mu_diff
        V =  U * V; % trace(AA')=sum(A(:).^2)
        trM = sum(V(:).^2); % trM = trace(Sigma_2^-1 *Sigma1)
        kl = 0.5*(trM + term - conddim) - logdetR;
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
        str = [str, sprintf('covariance (%d-by-%d): %s\n', size(theta.sigma,1), size(theta.sigma,2), mat2str(theta.sigma, 4))];
        
        if nargout == 0
            str = [sprintf('%s distribution parameters:\n', D.name()), str];
        end
    end

%%

    D = mxe_addsharedfields(D);
end
