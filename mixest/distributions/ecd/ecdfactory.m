%% |ecdfactory|
% Construct a zero-mean elliptically contoured distribution structure
%
% *Syntax*
%
%   D = ecdfactory(datadim, radialD)
%
% *Description*
%
% |D = ecdfactory(datadim, radialD)| returns a structure representing a
% |datadim|-dimensional ecd distribution with square radial pdf |radialD|.
%
% *Distribution Parameters*
%
% * *|sigma|* (|datadim-by-datadim|) The scatter matrix.
% * *|radialD|* (parameter structure) Contains the parameters of the square
% radial pdf
%
% *Probability Density Function*
%
% The distribution has the following density:
% 
% $$ f(x;\mu,\Sigma)=
% |\Sigma|^{-1/2} \pi^{-n/2} \gamma{n/2} (x^T \Sigma^{-1} x)^{1-n/2}
%   p(x^T \Sigma^{-1} x) $$
%
% where $n$ is the data space dimensions, $\Sigma$ is the scatter matrix
% and $p$ is the squared radial distribution.
% 
% *Example*
%
%   % Construct a 2-dimensional elliptical-gamma distribution
%   D = ecdfactory(2, gammafactory());
%   % Build a parameter structure for it:
%   theta = struct('sigma', [3 0; 0 2], 'radialD', struct('a', 3, 'b', 1));
%   % Plot the PDF:
%   x = -5:0.2:5;
%   y = -5:0.2:5;
%   [X, Y] = meshgrid(x, y);
%   data = [X(:) Y(:)]';
%   f = D.pdf(theta, data);
%   surf(X, Y, reshape(f, size(X)));
%
% <<img/ecdfactory_01.png>>
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

function D = ecdfactory(datadim, radialD, fixing)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'ecd';

%%

    assert(datadim >= 1, 'datadim must be an integer larger than or equal to 1.');

%% 
% Flag to control the memory usage (resulting code will be slower)

    low_memory = true;

%%

    thetamask = struct;
    fixed_radialD = false;
    fixed_sigma = false;
    fixed_data = false;
    fixed_u = [];
    fixed_logdetR = [];
    
    % if the first argument is a struct with a 'fixing__' field, we should
    % make some parameters fixed
    if nargin>2 && isfield(fixing, 'fixing__')
        % fixing.thetamask is a parameter structure with non-existing or
        % NaN values for variable parameters and the desired fixed values
        % for fixed parameters.
        thetamask = fixing.thetamask;
        if isfield(thetamask, 'radialD') 
            fixed_radialD = true;
        end
        if isfield(thetamask, 'sigma') 
            %disp('fixing sigma');
            fixed_sigma = true;
        end
        if isfield(fixing, 'data') && fixed_sigma
            %disp('fixing data');
            fixed_data = true;
            store = struct;
            store = ll_intermediate_params1(thetamask, store);
            store = ll_intermediate_params2(fixing.data, store);
            fixed_u = store.u;
            fixed_logdetR = store.logdetR;
        end
    end    
    
%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.

    D.M = build_param_manifold();
    function M = build_param_manifold()
        elements = struct;
        if ~fixed_radialD
            elements.radialD = radialD.M;
        end
        if ~fixed_sigma
            elements.sigma = spdfactory(datadim);
        end
        M = productmanifold(elements);
    end

%% |dim|
% See <doc_distribution_common.html#3 distribution structure common members>.

    D.dim = @() datadim*(datadim+1)/2 + datadim + radialD.dim() - 1; % parameter space dimensions

%% |datadim|
% See <doc_distribution_common.html#4 distribution structure common members>.

    D.datadim = @() datadim; % data space dimensions
    
%% |radialD|
% Returning radial distribution

    D.radialD = @() radialD;

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

    function store = ll_intermediate_params2(data, store)
    % calculate intermediate parameters for ll #2
    % Note: store must contain fields created by the previous ll_intermediate_params functions
        
        if ~isfield(store, 'u') || ~isfield(store, 'Rinvdata') 
            data = mxe_readdata(data);
            data = data.data;
            % u = X' C^-1 X
            store.Rinvdata = (store.Rinv' * data);
            store.u = sum(store.Rinvdata.^2, 1);
        end
        
        
    end

    function store = ll_intermediate_params3(theta, store)
    % calculate intermediate parameters for ll #3
    % Note: store must contain fields created by the previous ll_intermediate_params functions
        
        if ~isfield(store, 'radialDllvec')
            if ~isfield(store,'radialD')
                store.radialD = struct;
            end
            [store.radialDllvec, store.radialD] = radialD.llvec(...
                theta.radialD, store.u, store.radialD);
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
        
        if fixed_sigma && ~fixed_data
            warning('Sigma is fixed but data was not given during fixation');
        end
        
        theta = fullparam(theta);
        
        if ~fixed_data
            store = ll_intermediate_params1(theta, store);
            store = ll_intermediate_params2(data, store);
            logdetR = store.logdetR;
        else
            data = mxe_readdata(data, false);
            store.u = fixed_u(data.index);
            logdetR = fixed_logdetR;
        end
        
        store = ll_intermediate_params3(theta, store);
                    
        llvec = - (0.5*datadim)*log(pi) - logdetR + gammaln(0.5*datadim) ...
            + store.radialDllvec + (1-0.5*datadim)*log(store.u);
        
        if ~isempty(weight)
            llvec = weight .* llvec;
        end
        
        if low_memory
            store = rmfield(store, 'Rinvdata');
        end
        
    end

%% |llgrad|
% See <doc_distribution_common.html#7 distribution structure common members>.

    D.llgrad = @llgrad;
    function [dll, store] = llgrad(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        theta = fullparam(theta);
        
        weight = mxe_readweight(data);
        
        if fixed_sigma && ~fixed_data
            warning('Sigma is fixed but data was not given during fixation');
        end
        
        if ~fixed_data
            store = ll_intermediate_params1(theta, store);
            store = ll_intermediate_params2(data, store);
        else
            data = mxe_readdata(data, false);
            store.u = fixed_u(data.index);
        end
        
        store = ll_intermediate_params3(theta, store);

        u = store.u;
        
        if ~fixed_radialD
            [dll.radialD, store.radialD] = radialD.llgrad(theta.radialD, ...
                struct('data',u, 'weight', weight), store.radialD);
        end
        
        if ~fixed_sigma
            if ~isfield(store, 'radialDpart')
                [Dgraddata, store.radialD] = radialD.llgraddata(theta.radialD, u, store.radialD);
                store.radialDpart = Dgraddata + (1-0.5*datadim)./u;
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
        
        if low_memory && false
            store = rmfield(store, 'Rinvdata');
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
        
        store = ll_intermediate_params1(theta, store);
        
        store = ll_intermediate_params2(data, store);
        store = ll_intermediate_params3(theta, store);
        
        if ~isfield(store, 'radialDpart')
            [Dgraddata, store.radialD] = radialD.llgraddata(theta.radialD, store.u, store.radialD);
            store.radialDpart = Dgraddata + (1-0.5*datadim) ./ store.u;
        end
        
        dld =  2 * store.Rinv * store.Rinvdata;
        
        if ~isempty(weight)
            dld = bsxfun(@times, dld, weight.* store.radialDpart);
        else
            dld = bsxfun(@times, dld, store.radialDpart);
        end
        
        if low_memory && false
            store = rmfield(store, 'Rinvdata');
        end
        
    end

%% |gaussianize|
% 
    D.gaussianize = @gaussianize;
    function y = gaussianize(theta, data)
        store = struct;
        store = ll_intermediate_params1(theta, store);
        store = ll_intermediate_params2(data, store);
        % change the radial distribution to chi-squared
        u = chi2inv(radialD.cdf(theta.radialD, store.u), datadim); %TODO implement chi2inv
        y = bsxfun(@times, sqrt(u./store.u), store.Rinvdata);
    end

%% |pdf|
% See <doc_distribution_common.html#10 distribution structure common members>.

%% |sample|
% See <doc_distribution_common.html#11 distribution structure common members>.

    D.sample = @sample;
    function data = sample(theta, n)
        
        if nargin < 2, n = 1; end
        theta = fullparam(theta);
        data = randn(datadim, n);
        % sampling from squared radial component
        r2 = radialD.sample(theta.radialD, n); 
        % projective normal data into the unit sphere & multiply by radial
        data = bsxfun(@times, data , sqrt( r2 ./ sum(data.^2,1)) );
        L = chol(theta.sigma,'lower');
        data = L * data;

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
        if ~fixed_sigma
            if isempty(weight)
                sigma = 1/N * (data * data.');
            else
                data2 = bsxfun(@times, sqrt(weight), data);
                sigma = 1/sum(weight) * (data2 * data2.');
            end
            theta.sigma = sigma;
        else
            sigma = thetamask.sigma;
        end
        R = chol(sigma);
        Rinv = R \ eye(datadim);
        Rinvdata = (Rinv' * data);
        u = sum(Rinvdata.^2, 1);
        if ~fixed_radialD
            theta.radialD = radialD.init(u);
        end
    end

%%

% %% |estimatedefault|
     D.estimatedefault = @estimatedefault;
     function [varargout] = estimatedefault(varargin)
         data = varargin{1};
         if fixed_sigma && ~fixed_data
            warning('Sigma is fixed but data was not given during fixation');
            store = struct;
            store = ll_intermediate_params1(thetamask, store);
            store = ll_intermediate_params2(data, store);
            fixed_u = store.u;
         end
         if fixed_sigma
             if nargin > 1
                 options = varargin{2};
                 theta0 = options.theta0;
                 if ~isempty(theta0)
                     options.theta0 = theta0.radialD;
                     varargin{2} = options;
                 end
             end
             % using default estimator of the radial component
             data = mxe_readdata(data, false);
             data.data = fixed_u;
             varargin{1} = data;
             [varargout{1:nargout}] = radialD.estimatedefault(varargin{:});
             theta = varargout{1};
             varargout{1} = struct;
             varargout{1}.radialD = theta;
         elseif fixed_radialD && strcmp (radialD.name(), 'gamma')
             [varargout{1:nargout}] = ecd_eg_estimatedefault(D, ...
                 thetamask.radialD, varargin{:}); 
         end
         if ~fixed_sigma && ~fixed_radialD 
             [varargout{1:nargout}] = ecd_estimatedefault(D, varargin{:});
         end
     end

%% |penalizerparam|
% See <doc_distribution_common.html#15 distribution structure common members>.
%
% *Penalizer Info*
%
% The default penalizer for this distribution is the Inverse-Wishart
% distribution for covariance and default prior for squared-radial.
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

    D.penalizerparam = @penalizerparam;
    function penalizer_theta = penalizerparam(data)
        % Using Fraley&Raftery (2007) method for computing parameters
        data = mxe_readdata(data);
        n = data.size;
        data = data.data;
        
        penalizer_theta.nu = 2;
        sigma = (data * data.') / n;
        % Check if it is multiplication of identity
        if isequal(sigma, sigma(1,1) * eye(datadim))
            penalizer_theta.invLambda = 1/sigma(1,1);
        else
            penalizer_theta.invLambda = inv(sigma);
        end
        % Computing squared-radial distribution
        R = chol(sigma);
        Rinv = R \ eye(datadim);
        Rinvdata = (Rinv' * data);
        u = sum(Rinvdata.^2, 1);
        penalizer_theta.radialD = radialD.penalizerparam(u);
    end

%% |penalizercost|
% See <doc_distribution_common.html#16 distribution structure common members>.

    D.penalizercost = @penalizercost;
    function [costP, store] = penalizercost(theta, penalizer_theta, store)
        
        if nargin < 3
            store = struct;
        end
        theta = fullparam(theta);
        if ~isfield(store,'radialD')
            store.radialD = struct;
        end
        
        store = ll_intermediate_params1(theta, store);
        logdetR = store.logdetR;
        
        if ~isfield(store, 'Sinv')
            store.Sinv = store.Rinv * store.Rinv';
        end
            
        Sinv = store.Sinv;

        const2 = penalizer_theta.nu + datadim + 1;
        
        [costP, store.radialD] = radialD.penalizercost(theta.radialD, penalizer_theta.radialD, store.radialD);
        
        % If data to penalizeparam is whitened then invLambda is scalar
        if isscalar(penalizer_theta.invLambda)
            costP = costP - const2 * logdetR - 0.5 * penalizer_theta.invLambda * trace(Sinv);
        else
            costP = costP - const2 * logdetR - 0.5 * (penalizer_theta.invLambda(:).' * Sinv(:));
        end
    end

%% |penalizergrad|
% See <doc_distribution_common.html#17 distribution structure common members>.

    D.penalizergrad = @penalizergrad;
    function [gradP, store] = penalizergrad(theta, penalizer_theta, store)
        
        if nargin < 3
            store = struct;
        end
        theta = fullparam(theta);
        if ~isfield(store,'radialD')
            store.radialD = struct;
        end
        
        store = ll_intermediate_params1(theta, store);
        Sinv = store.Sinv;

        const2 = penalizer_theta.nu + datadim + 1;
        
        const2 = - 0.5 * const2;
        
        gradP.radialD = radialD.penalizergrad(theta.radialD, penalizer_theta.radialD, store.radialD);
        
        if isscalar(penalizer_theta.invLambda)
            gradP.sigma = const2 * Sinv +  ...
                (0.5 * penalizer_theta.invLambda) * (Sinv * Sinv.');
        else
            gradP.sigma = const2 * Sinv +  ...
                0.5 * Sinv * penalizer_theta.invLambda * Sinv;
        end

    end

%% |sumparam|
% See <doc_distribution_common.html#18 distribution structure common members>.

    D.sumparam = @sumparam;
    function theta = sumparam(theta1, theta2)
        if ~fixed_radialD
            theta.radialD = radialD.sumparam(theta1.radialD , theta2.radialD);
        end
        if ~fixed_sigma
            theta.sigma = theta1.sigma + theta2.sigma;
        end
    end

%% |scaleparam|
% See <doc_distribution_common.html#19 distribution structure common members>.

    D.scaleparam = @scaleparam;
    function theta = scaleparam(scalar, theta)
        if ~fixed_radialD
            theta.radialD = radialD.scaleparam(scalar , theta.radialD);
        end
        if ~fixed_sigma
            theta.sigma = scalar * theta.sigma;
        end
    end

%% |sumgrad|
% See <doc_distribution_common.html#20 distribution structure common members>.

%% |scalegrad|
% See <doc_distribution_common.html#21 distribution structure common members>.

%% |entropy|
% See <doc_distribution_common.html#22 distribution structure common members>.

    D.entropy = @entropy;
    function h = entropy(theta)
        theta = fullparam(theta);
        logdetR = sum(log( diag( chol(theta.sigma) ) )); 

        function y = intfunc(v)
            lvec = radialD.llvec(theta.radialD, v);
            y = ((1 - 0.5 * datadim) * log(v) + lvec) .* exp(lvec);
        end
        
        h = logdetR  + 0.5*datadim*log(pi) - gammaln(0.5 * datadim)- ...
          quadgk(@(x)intfunc(x), 0, Inf);          
    end

%% |kl|
% See <doc_distribution_common.html#23 distribution structure common members>.

    D.kl = @kl;
    function kl = kl(theta1, theta2) %#ok<STOUT,INUSD>
        error('KL-divergence not implemented')
    end

%% |AICc|
% See <doc_distribution_common.html#24 distribution structure common members>.

%% |BIC|
% See <doc_distribution_common.html#25 distribution structure common members>.

%% |display|
% See <doc_distribution_common.html#26 distribution structure common members>.

    D.display = @display;
    function str = display(theta)
        theta = fullparam(theta);
        str = '';
        str = [str, sprintf('sigma (%d-by-%d): %s\n', size(theta.sigma,1), size(theta.sigma,2), mat2str(theta.sigma, 4))];
        str = [str, sprintf('Radial Distribution: %s\n', radialD.name())];
        str = [str, radialD.display(theta.radialD)];
        
        if nargout == 0
            str = [sprintf('%s distribution parameters:\n', D.name()), str];
        end
    end

%% |fixate|
% Make some parameters fixed
%
% *Syntax*
%
%   newD = D.fixate(paramName, paramValue, ...)
%
% *Description*
%
% |newD = D.fixate(paramName, paramValue, ...)| returns the new gamma
% distribution |newD| where the specified parameters are fixed so that they
% won't change in estimation. |paramName| is the name of the parameter ('a'
% or 'b') and |paramValue| is its desired fixed value.
%
% *Example*
%
%   D = gammafactory();
%   D = D.fixate('a', 2);
% 


    D.fixate = @fixate;
    function newD = fixate(varargin)
    % syntax: newD = fixate('a', 1, 'b', 2, ...)
        
        % build theta mask
        Athetamask = struct(varargin{:});
        
        % build fixing struct
        Afixing = struct('fixing__', true);
        if isfield(Athetamask, 'data');
            Afixing.data = Athetamask.data;
            Athetamask = rmfield(Athetamask, 'data');
        end
        Afixing.thetamask = Athetamask;
        
        % construct the new distribution
        newD = ecdfactory(datadim, radialD, Afixing);
    end

%% |unfix|
% Undo |fixate|
%
% *Syntax*
%
%   newD = D.unfix(paramName, ...)
%
% *Description*
%
% |newD = D.unfix(paramName, ...)| where |D| is a gamma distribution with
% fixed parameters, returns the new gamma distribution |newD| where the
% specified parameters are unfixed. |paramName| is the name of the
% parameter ('a' or 'b').
%
% *Example*
%
%   D = gammafactory();
%   D = D.fixate('a', 2, 'b', 1);
%   D = D.unfix('a');
% 

    D.unfix = @unfix;
    function newD = unfix(varargin)
    % syntax: newD = unfix('a', 'b', ...)
    
        % build theta mask
        Athetamask = rmfield(thetamask, varargin);
        
        % build fixing struct
        Afixing = struct('fixing__', true);
        Afixing.thetamask = Athetamask;
        
        % construct the new distribution
        newD = ecd3factory(datadim, radialD, Afixing);
    end

%% |fullparam|
% Get the full parameter structure, given the variable parameter structure
%
% *Syntax*
%
%   fulltheta = D.fullparam(theta)
%
% *Description*
%
% |fulltheta = D.fullparam(theta)| where |D| is a gamma distribution with
% fixed parameters, returns the full parameters of the gamma distribution
% including the fixed parameters. |theta| is a parameter structure for
% the distribution |D| (i.e. excluding the fixed parameters).
%
% *Example*
%
%   D = gammafactory();
%   D = D.fixate('a', 2);
%   theta = struct('b', 1);
%   fulltheta = D.fullparam(theta)
%
%  fulltheta = 
% 
%      b: 1
%      a: 2
%

    D.fullparam = @fullparam;
    function fulltheta = fullparam(theta)
        
        fulltheta = theta;
        if fixed_sigma
            fulltheta.sigma = thetamask.sigma;
        end
        if fixed_radialD
            fulltheta.radialD = thetamask.radialD;
        end
    end

%%

    D = mxe_addsharedfields(D);
end