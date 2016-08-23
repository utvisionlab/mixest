%% |agfactory|
% Construct an angular-Gaussian distribution structure
%
% *Syntax*
%
%   D = agfactory(datadim)
%
% *Description*
%
% |D = agfactory(datadim)| returns a structure representing a
% |datadim|-dimensional angular-Gaussian distribution.
%
% *Distribution Parameters*
%
% * *|sigma|* (|datadim-by-datadim| matrix)
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Contributors:
%  Reshad Hosseini
%  Mohamadreza Mash'al
%  Poorya Habibzadeh
%
% Change log: 
%

function D = agfactory(datadim)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'ag';

%%

    assert(datadim >= 1, 'datadim must be an integer larger than or equal to 1.');

%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.

    sigma = spdfactory(datadim);
    D.M = mxe_productmanifold(struct('sigma', sigma));

%% |dim|
% See <doc_distribution_common.html#3 distribution structure common members>.

    D.dim = @() datadim*(datadim+1)/2 ; % parameter space dimensions

%% |datadim|
% See <doc_distribution_common.html#4 distribution structure common members>.

    D.datadim = @() datadim; % data space dimensions

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
                    theta.sigma = theta.sigma + t * eye(datadim);
                    [R, p] = chol(theta.sigma); %#ok<ASGLU>
                    t = t * 100;
                end
                R = chol(theta.sigma);
            end
            store.Rinv = R \ eye(datadim); % Faster version of Rinv = inv_triu(R);
            store.logdetR = sum(log(diag(R)));
        end

    end

    function store = ll_intermediate_params2(theta, data, store)
    % calculate intermediate parameters for ll #2
    % Note: store must contain fields created by the previous ll_intermediate_params functions
        
        if ~isfield(store, 'u') || ~isfield(store, 'Rinvdata') || ~isfield(store, 'Dllvec')
            data = mxe_readdata(data);
            data = data.data;
            % u = X' C^-1 X
            store.Rinvdata = (store.Rinv' * data);
            store.u = sum(store.Rinvdata.^2, 1);
            %if ~isfield(store,'radialD')
              %  store.radialD = struct;
            %end
            %[store.radialDllvec, store.radialD] = radialD.llvec(theta.radialD, store.u, store.radialD);
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
        
        weight = mxe_readweight(data);
        
        store = ll_intermediate_params1(theta, store);
        logdetR = store.logdetR;

        store = ll_intermediate_params2(theta, data, store);
                    
        llvec = - (0.5*datadim)*log(pi) - log(2) - logdetR + gammaln(0.5*datadim) ...
            -0.5*datadim*log(store.u);
        
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
        
        weight = mxe_readweight(data);
        
        store = ll_intermediate_params1(theta, store);
        
        store = ll_intermediate_params2(theta, data, store);

        u = store.u;
        
         if ~isfield(store, 'radialDpart')
             store.radialDpart =  -0.5*datadim./u;
         end

        if ~isempty(weight)
            w = store.radialDpart .* weight ;
        else
            w =  store.radialDpart ;
        end
        
        if ~isempty(weight)
            sumW = sum(weight, 2);
        else
            sumW = length(u);
        end
        
        Rinvdata = store.Rinvdata;
        Rinv = store.Rinv;
        Rinvdata2 = bsxfun(@times, w , Rinvdata);
        dsigma =  - 0.5 * sumW * eye(datadim) - (Rinvdata2 * Rinvdata.');
        dll.sigma =  Rinv * dsigma * Rinv.';
        
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
        
        store = ll_intermediate_params2(theta, data, store);
        
         u = store.u;
         
        if ~isfield(store, 'radialDpart')
            store.radialDpart =  -0.5*datadim./u;
        end
        
        dld =  2 * store.Rinv * store.Rinvdata;
        
        if ~isempty(weight)
            dld = bsxfun(@times, dld, weight.* store.radialDpart);
        else
            dld = bsxfun(@times, dld, store.radialDpart);
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
        L = chol(theta.sigma,'lower');
        data = L * data;
        data = bsxfun(@rdivide, data, sqrt(sum(data.^2,1)));
    end

%% |randparam|
% See <doc_distribution_common.html#12 distribution structure common members>.

%% |init|
% See <doc_distribution_common.html#13 distribution structure common members>.

    D.init = @init;
    function theta = init(data, varargin)
        data = mxe_readdata(data);
        weight = data.weight;
        N = data.size;
        data = data.data;

        if isempty(weight)
            sigma = 1/N * (data * data.');
        else
            data2 = bsxfun(@times, sqrt(weight), data);
            N = sum(weight);
            sigma = 1/N * (data2 * data2.');
        end
        if N < datadim
            disp('number of data is small');
            sigma = sigma + mean(diag(sigma))/10*eye(datadim);
        end
        theta.sigma = sigma;
    end

%%

%% |estimatedefault|
    D.estimatedefault = @estimatedefault;
    function [varargout] = estimatedefault(varargin)
        [varargout{1:nargout}] = ag_estimatedefault(D, varargin{:}); 
    end

%% |penalizerparam|
% See <doc_distribution_common.html#15 distribution structure common members>.
%
% *Penalizer Info*
%
% The default penalizer leads to regularized Tyler estimator
%
% The form of penalizer:
%
% $$ f(\Sigma) = 
% -0.5*log(|\Sigma|) - 0.5*d* log(trace(\Sigma^{-1} invLambda))) $$
%
% where
%
% * *|invLambda|* (|datadim-by-datadim| matrix) : The inverse scale matrix.
%
 
    D.penalizerparam = @penalizerparam;
    function penalizer_theta = penalizerparam(data)
        data = mxe_readdata(data);
        N = data.size;
        data = data.data;
        
        if N < datadim
            disp('number of data is small in penalizer');
            sigmat = D.init(data);
        else
            sigmat = D.estimate(data);
        end
        
        sigmat = sigmat.sigma;
        if isequal(sigmat, sigmat(1,1) * eye(datadim))
            penalizer_theta.invLambda = sigmat(1,1);
        else
            penalizer_theta.invLambda = sigmat;
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
        
        if ~isfield(store, 'Sinv')
            Rinv = store.Rinv;
            store.Sinv = Rinv * Rinv';
        end
        
        Sinv = store.Sinv;
        
        costP = - logdetR;
        % If data to penalizeparam is whitened then invLambda is scalar
        if isscalar(penalizer_theta.invLambda)
            costP = costP - 0.5 * datadim * penalizer_theta.invLambda * log(trace(Sinv));
        else
            costP = costP - 0.5 * datadim* log(penalizer_theta.invLambda(:).' * Sinv(:));
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
        
        if ~isfield(store, 'Sinv')
            Rinv = store.Rinv;
            store.Sinv = Rinv * Rinv';
        end
        
        Sinv = store.Sinv;
        
        if isscalar(penalizer_theta.invLambda)
            gradP.sigma = -0.5 * Sinv +  0.5 * datadim * ...
                penalizer_theta.invLambda * (Sinv * Sinv.');
        else
            slambda = (penalizer_theta.invLambda(:).' * Sinv(:));
            gradP.sigma = -0.5 * Sinv +  0.5 * datadim * ...
                1./slambda * Sinv * penalizer_theta.invLambda * Sinv;
        end
 
    end


%% |sumparam|
% See <doc_distribution_common.html#18 distribution structure common members>.

    D.sumparam = @sumparam;
    function theta = sumparam(theta1, theta2)
        theta.sigma = theta1.sigma + theta2.sigma;
    end

%% |scaleparam|
% See <doc_distribution_common.html#19 distribution structure common members>.

    D.scaleparam = @scaleparam;
    function theta = scaleparam(scalar, theta)
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
        mat = eig(theta.sigma);
        logdet = sum(log(mat));
        ELog = Integral_Method(mat);
        h = 0.5 * logdet - (datadim/2) * (ELog - psi(datadim/2) - log(2)) - ...
            gammaln(datadim/2) + (datadim/2) * log(pi) + log(2);
    end

%% |kl|
% See <doc_distribution_common.html#23 distribution structure common members>.

    D.kl = @kl;
    function kl = kl(theta1, theta2) %#ok<STOUT,INUSD>
        mat = eig(theta1.sigma, theta2.sigma);
        logdet = sum(log(mat));
        ELog = Integral_Method(mat);
        kl = -0.5 * logdet + (datadim/2) * (ELog - psi(datadim/2) - log(2));
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
        str = [str, sprintf('sigma (%d-by-%d): %s\n', size(theta.sigma,1), size(theta.sigma,2), mat2str(theta.sigma, 4))];
        
        if nargout == 0
            str = [sprintf('%s distribution parameters:\n', D.name()), str];
        end
    end

%%

    D = mxe_addsharedfields(D);
end