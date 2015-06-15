%% |cmvn2factory|
% Construct a reparameterized conditional multi-variate normal distribution
% structure
%
% *Syntax*
%
%   D = cmvn2factory(datadim, conddim)
%
% *Description*
%
% |D = cmvnfactory(datadim, conddim)| returns a structure representing a
% conditional normal distribution with a special reparameterization.
% |datadim| is the data dimensions and |conddim| is the condition
% dimensions.
%
% *See also* <cmvnfactory.html cmvnfactory>
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

function D = cmvn2factory(datadim, conddim)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'cmvn2';

%%

    assert(datadim >= 1, 'datadim must be an integer larger than or equal to 1.');

%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.

    %W = euclideanfactory(conddim, datadim+1);
    sigmaT = spdfactory(conddim + datadim + 1);
    D.M = productmanifold(struct('sigmat', sigmaT));

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

% Joint and marginal manifold
%    
    Djoint = mvn2factory(conddim + datadim);
    Dmarg = mvn2factory(datadim);
    
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
            theta.sigmat = zeros(datadim+1);
            theta.sigmat(1:datadim+1,datadim+2:end) = theta.W;
            theta.sigmat(datadim+2:end,1:datadim+1) = theta.W.';
            theta.sigmat(datadim+2:end,datadim+2:end) = theta.sigma + theta.W*theta.W.';
            theta.sigmat(1:datadim+1,1:datadim+1) = eye(datadim+1);
            theta = rmfield(theta,'mu');
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
        if isfield(theta,'sigmat');
            sigma = theta.sigmat;
            sigma1 = sigma(datadim+2:end,datadim+2:end );
            sigma2 = sigma(1:datadim+1, 1:datadim+1);
            sigma12 = sigma(datadim+2:end, 1:datadim+1);
            theta.W = sigma12 / sigma2;
            theta.sigma = sigma1 - theta.W * sigma12.';
            theta = rmfield(theta,'sigmat');
        end
    end

%% |ll|
% See <doc_distribution_common.html#5 distribution structure common members>.

    D.ll = @ll;
    function [ll, store] = ll(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        if ~isfield(store,'joint');
            store.joint = struct;
        end
        
        if ~isfield(store,'marg');
            store.marg = struct;
        end
        
        data = mxe_readdata(data);
        indjoint = [1:datadim, datadim+2:datadim+conddim+1, datadim+1];
        theta.sigmat = theta.sigmat(indjoint, indjoint);
        [lljoint, store.joint] = Djoint.ll(theta, data, store.joint);
        
        indmarg = [1:datadim datadim+conddim+1];
        theta.sigmat = theta.sigmat(indmarg, indmarg);
        data.data = data.data(1:datadim ,:);
        [llmarg, store.marg] = Dmarg.ll(theta, data, store.marg);
         
        ll = lljoint - llmarg;
    end

%% |llvec|
% See <doc_distribution_common.html#6 distribution structure common members>.

    D.llvec = @llvec;
    function [llvec, store] = llvec(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        if ~isfield(store,'joint');
            store.joint = struct;
        end
        
        if ~isfield(store,'marg');
            store.marg = struct;
        end
        
        data = mxe_readdata(data);
        indjoint = [1:datadim, datadim+2:datadim+conddim+1, datadim+1];
        theta.sigmat = theta.sigmat(indjoint, indjoint);
        [llvjoint, store.joint] = Djoint.llvec(theta, data, store.joint);
        
        indmarg = [1:datadim datadim+conddim+1];
        theta.sigmat = theta.sigmat(indmarg, indmarg);
        data.data = data.data(1:datadim ,:);
        [llvmarg, store.marg] = Dmarg.llvec(theta, data, store.marg);
         
        llvec = llvjoint - llvmarg;
    end

%% |llgrad|
% See <doc_distribution_common.html#7 distribution structure common members>.

    D.llgrad = @llgrad;
    function [dll, store] = llgrad(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        if ~isfield(store,'joint');
            store.joint = struct;
        end
        
        if ~isfield(store,'marg');
            store.marg = struct;
        end
        
        data = mxe_readdata(data);
        indjoint = [1:datadim, datadim+2:datadim+conddim+1, datadim+1];
        theta.sigmat = theta.sigmat(indjoint, indjoint);
        [dlljoint, store.joint] = Djoint.llgrad(theta, data, store.joint);
        
        indmarg = [1:datadim datadim+conddim+1];
        theta.sigmat = theta.sigmat(indmarg, indmarg);
        data.data = data.data(1:datadim ,:);
        [dllmarg, store.marg] = Dmarg.llgrad(theta, data, store.marg);
         
        indjoint = [1:datadim, datadim+conddim+1, datadim+1:datadim+conddim];
        indmarg = 1:datadim+1;
        dll.sigmat = dlljoint.sigmat(indjoint, indjoint);
        dll.sigmat(indmarg, indmarg) = dll.sigmat(indmarg, indmarg) - dllmarg.sigmat;
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
        theta = new2main(theta);
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
        theta = new2main(theta);
        if isscalar(datamargin)
            n = datamargin;
            datamargin = randn(datadim, n);
        end
        n = size(datamargin,2);
        data = randn(conddim, n);
        L = chol(theta.sigmat,'lower');
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
%   % create a cmvn2 distribution
%   D = cmvn2factory(1, 1);
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
        sumW = n;
        if ~isempty(weight)
            data = bsxfun(@times, sqrt(weight), data);
            sumW = sum(weight);
        end
        sigmat = (data * data.')/sumW;
        if options.penalize
            denum = (n + penalizer_est_theta.nu);
            indjoint = [1:datadim, datadim+2:datadim+conddim+1, datadim+1];
            if isscalar(penalizer_est_theta.invLambda)
                 theta.sigmat = n/ denum * sigma + ...
                     1/denum * diag( repmat(penalizer_est_theta.invLambda, [1 size(sigmat,1)]));
            else
                theta.sigmat = n/ denum * sigmat + ...
                    1/denum * penalizer_est_theta.invLambda(indjoint, indjoint);
            end
        else 
            theta.sigmat = sigmat;
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
        penalizer_theta = Djoint.penalizerparam(data);
    end

%% |penalizercost|
% See <doc_distribution_common.html#16 distribution structure common members>.

    D.penalizercost = @penalizercost;
    function [costP, store] = penalizercost(theta, penalizer_theta, store)
        
        if nargin < 3
            store = struct;
        end
        
        if ~isfield(store,'joint');
            store.joint = struct;
        end
        
        [costP, store.joint] = Djoint.penalizercost(theta, penalizer_theta, store.joint);

    end

%% |penalizergrad|
% See <doc_distribution_common.html#17 distribution structure common members>.

    D.penalizergrad = @penalizergrad;
    function [gradP, store] = penalizergrad(theta, penalizer_theta, store)
        
        if nargin < 3
            store = struct;
        end
        
         if ~isfield(store,'joint');
            store.joint = struct;
        end
        
        [gradP, store.joint] = Djoint.penalizergrad(theta, penalizer_theta, store.joint);

    end

%% |sumparam|
% See <doc_distribution_common.html#18 distribution structure common members>.

    D.sumparam = @sumparam;
    function theta = sumparam(theta1, theta2)
        theta.sigmat = theta1.sigmat + theta2.sigmat;
    end

%% |scaleparam|
% See <doc_distribution_common.html#19 distribution structure common members>.

    D.scaleparam = @scaleparam;
    function theta = scaleparam(scalar, theta)
        theta.sigmat = scalar * theta.sigmat;
    end

%% |sumgrad|
% See <doc_distribution_common.html#20 distribution structure common members>.

%% |scalegrad|
% See <doc_distribution_common.html#21 distribution structure common members>.

%% |entropy|
% See <doc_distribution_common.html#22 distribution structure common members>.

    D.entropy = @entropy;
    function h = entropy(theta)
        theta = new2main(theta);
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
