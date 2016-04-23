%% |vmffactory|
% Construct a von Mises-Fisher distribution structure
%
% *Syntax*
%
%   D = vmffactory(datadim)
%
% *Description*
%
% |D = vmffactory(datadim)| returns a structure representing a
% |datadim|-dimensional VMF distribution on the unit sphere.
%
% *Distribution Parameters*
%
% * *|mu|* (|datadim-by-1| vector) : The mean vector.
% * *|kappa|* (scaler) : The kappa parameter.
%
% *Probability Density Function*
%
% The distribution has the following density:
% 
% $$ f(x;\mu,\kappa)= c_n(\kappa) exp(\kappa \mu^{T} x) $$
%
% where $n$ is the data space dimensions, $\mu$ is the mean vector and
% $c_n(\kappa)$ is the normalization constant given by:
%
% $$ c_n(\kappa) = \frac{\kappa^{n/2-1}}{(2\pi)^{n/2} I_{n/2-1}(\kappa)} $$
%
% where $I_r$ is the modified Bessel function of the first kind and order
% $r$.
% 
% *Example*
%
%   % Construct several vmf distribution on a sphere:
%   D = vmffactory(3);
%   % Build three parameter structures and visualize samples
%   colors = {[ 0 0 1], [0 1 0], [1 0 0]};
%   mus = {[-0.57; 0.57; 0.57], [0; -1; 0], [1; 0; 0]};
%   for k = 1 : 3
%       theta.kappa = 7^(k-1)*2;
%       theta.mu = mus{k}; 
%       data = D.sample(theta, 5000);
%       h = plot3(data(1,:), data(2,:), data(3,:), '.');
%       set(h, 'MarkerSize', 4, 'Color', colors{k}, 'Marker', '.');
%       h = line([0 mus{k}(1)], [0 mus{k}(2)], [0 mus{k}(3)]);
%       set(h, 'LineWidth', 2, 'Color', colors{k});
%       hold on
%   end
%
% <<img/vmffactory_01.png>>
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

function D = vmffactory(datadim)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'vmf';

%%

    assert(datadim >= 2, 'datadim must be an integer larger than or equal to 2.');

%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.

    muM = spherefactory(datadim);
    kappaM = positivefactory();
    D.M = productmanifold(struct('mu', muM, 'kappa', kappaM));

%% |dim|
% See <doc_distribution_common.html#3 distribution structure common members>.

    D.dim = @() datadim; % parameter space dimensions

%% |datadim|
% See <doc_distribution_common.html#4 distribution structure common members>.

    D.datadim = @() datadim; % data space dimensions3
    
%%

    function [I, store] = bessln(theta, store)
    % besseliln(nu,z) is the log modified bessel function of the first kind
        
        if nargin < 3
            store = struct;
        end
        
        nu = datadim/2 - 1;
        x = theta.kappa;
        if false
            % More accurate but needs integration
            if nu == 0
                % First method applicable for nu = 0
                fun = @(t) exp(x*(cos(t)-1));
                store.Iint = quadgk(fun, 0, pi, 'RelTol', 1e-6, 'AbsTol', 0);
                I = x - log(pi) + log(store.Iint);
            elseif nu == 0.5
                % Second method applicable for nu = 0.5
                store.Iint = exp(-2*x);
                I = - 0.5 * log (x*2) - 0.5 * log(pi) + x + log(1 - store.Iint);
            else
                % Third method applicable for nu larger than 0.5
                v = nu - 0.5;
                a_x = -x ./(v + sqrt(v^2 + x.^2));
                h_x = a_x .* x - v * log(1 - a_x.^2);
                store.h_x = h_x;
                fun = @(t) exp(h_x - x * t + v * log(1 - t.^2) );
                store.Iint = quadgk(fun, -1 , a_x, 'RelTol', 1e-6, 'AbsTol', 0) + ...
                    quadgk(fun, a_x , 1, 'RelTol', 1e-6, 'AbsTol', 0);
                I = nu * log (x/2) - 0.5 * log(pi) - gammaln(v+1) - h_x + log(store.Iint);
            end
        else
            % Faster Code based on Abramowitz and Stegun
            if nu == 0
                t = x / 3.75;
                if t < 1
                    I = 1 + 3.5156229 * t.^2 + 3.0899424 * t.^4 + ...
                        1.2067492 * t.^6 + 0.2659732 * t.^8 + 0.0360768 ...
                        * t.^10 + 0.0045813 * t.^12;
                    store.Iint = I;
                    I = log(I);
                else
                    t = 1 / t;
                    I = 0.39894228 + 0.01328592 * t + 0.00225319 * t.^2 ...
                       - 0.00157565* t.^3 + 0.00916281 * t.^4 -  ...
                       0.02057706 * t.^5 + 0.02635537 * t.^6 - 0.01647633 ...
                       * t.^7 + 0.00392377 * t.^8;
                   store.Iint = I;
                   I = log(I) + x - 0.5 * log(x);
                end
            else
                frac = x/nu;
                square = 1 + frac.^2;
                root = sqrt(square);
                eta = root + log(frac) - log(1+root);
                I = - log(sqrt(2*pi*nu)) + nu*eta - 0.25*log(square);
            end
        end
        if isinf(I) || isnan(I)
            store.inf = true;
            warning('Inf or NaN bessln')
        else 
            store.inf = false;
        end
        
    end

    function [dI, store] = besslngrad(theta, store)
    % besslngrad(nu,z) is the gradient of the bessln with respect to kappa
    
    nu = datadim/2 - 1;
    x = theta.kappa;
    if false
        % More accurate but needs integration
        if ~isfield(store,'Iint')
            [notused, store] = bessln(theta, store);
        end
        if nu == 0
            % First method applicable for nu = 0
            %t2 = (cos(t)-1);
            fun = @(t) exp(x.*(cos(t)-1)) .* (cos(t)-1);
            Iint2 = quadgk(fun, 0, pi, 'RelTol', 1e-6, 'AbsTol', 0);
            dI = 1 + Iint2 / store.Iint;
        elseif nu == 0.5
            % Second method applicable for nu = 0.5
            dI = -0.5 / x + 1 + 2 * store.Iint / (1 - store.Iint);
        else
            % Third method applicable for nu larger than 0.5
            v = nu - 0.5;
            h_x = store.h_x;
            fun2 = @(t) t .* exp(h_x - x * t + v * log(1 - t.^2) );
            Iint2 = quadgk(fun2, -1 , 1,'RelTol', 1e-6, 'AbsTol', 0);
            dI = nu/x - Iint2/store.Iint;
        end
    else
        % Faster Code based on Abramowitz and Stegun
        if nu == 0
            if ~isfield(store,'Iint')
                [notused, store] = bessln(theta, store);
            end
            t = x / 3.75;
            if t < 1
                dI = 7.0312458 * t + 12.3597696 * t.^3 + ...
                    7.2404952 * t.^5 + 2.1277856 * t.^7 + 0.3607680 ...
                    * t.^9 + 0.0549756 * t.^11;
                dI = 1 / 3.75 ./ store.Iint .* dI ;
            else
                t = 1 / t;
                dI = 0.01328592 + 0.00450638 * t ...
                    - 0.00472695 * t.^2 + 0.03665124 * t.^3 -  ...
                    0.10288530 * t.^4 + 0.15813222 * t.^5 - 0.11533431 ...
                    * t.^6 + 0.03139016 * t.^7;
                dI = 1 - 1 / 3.75 ./store.Iint .* t.^2 .* dI - 0.5 / x;
            end
        else
            frac = x/nu;
            square = 1 + frac.^2;
            root = sqrt(square);
            dI = -0.5 * frac /square /nu + 1 /frac +  frac/(1+root);
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
        
        data = mxe_readdata(data);
        weight = data.weight;
        data = data.data;
        
        if ~isfield(store,'mux')
            store.mux = theta.mu.' * data;
        end
        
        mux = store.mux; 
        
        [I, store] = bessln(theta, store);
        
        const =  (datadim/2-1)*log(theta.kappa) - (datadim/2)*log(2*pi) - I;
        
        llvec = const + theta.kappa * mux;
        
        if ~isempty(weight)
            llvec = weight .* llvec;
        end
        if store.inf
            llvec(isinf(llvec) | isnan(llvec)) = -1e200;
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
        n = data.size;
        data = data.data;
        
        [dI, store] = besslngrad(theta, store);
        
        if ~isfield(store,'mux')
            store.mux = theta.mu.' * data;
        end
        
        if ~isempty(weight)
            sumW = sum(weight, 2);
            mux = store.mux .* weight;
            data = bsxfun(@times, weight, data);
        else
            sumW = n;
            mux = store.mux;
        end
       
        dll.mu = theta.kappa * sum(data, 2);
        
        dll.kappa = sumW * ((datadim/2-1)/theta.kappa - dI) + sum(mux) ;
        if isinf(dll.kappa) || isnan(dll.kappa)
            dll.kappa = 0;
        end

    end

%% |llgraddata|
% See <doc_distribution_common.html#8 distribution structure common members>.

    D.llgraddata = @llgraddata;
    function [dld, store] = llgraddata(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
 
        data = mxe_readdata(data);
        weight = data.weight;
        n = data.size;
        
        if ~isempty(weight)
            dld = bsxfun(@times, weight, theta.kappa * mu);
        else
            dld = repmat(theta.kappa * mu, [1 n]);
        end        
        
    end

%% |pdf|
% See <doc_distribution_common.html#10 distribution structure common members>.

%% |sample|
% See <doc_distribution_common.html#11 distribution structure common members>.

    D.sample = @sample;
    function data = sample(theta, n)
        % code was taken from vmf code by Suvrit and colleages
        data = vsamp(theta.mu, theta.kappa, n).';
    end

%% |randparam|
% See <doc_distribution_common.html#12 distribution structure common members>.

%% |init|
% See <doc_distribution_common.html#13 distribution structure common members>.

    D.init = @init;
    function theta = init(data, varargin)
        theta = estimatedefault(data);
        %theta.mu = muM.rand();
        %theta.kappa = 1 + randn(1)*1e-5;
    end

%% |estimatedefault|
% Default estimation function for von Mises-Fisher distribution. This
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
%   % create a VMF distribution
%   D = vmffactory(3);
%   % generate 1000 sample data points from the distribution
%   theta_sample = struct('mu', [0.7; 0.14; -0.7], 'kappa', 10);
%   data = D.sample(theta_sample, 1000);
%   % estimate distribution parameters to fit the data
%   theta = D.estimatedefault(data)
%

    D.estimatedefault = @estimatedefault;
    function [varargout] = estimatedefault(varargin)
        [varargout{1:nargout}] = vmf_estimatedefault(D, varargin{:});
    end

%% |penalizerparam|
% See <doc_distribution_common.html#15 distribution structure common members>.
%
% *Penalizer Info*
%
% The default penalizer for this distribution is the vMF 
% distribution for mean and Log-Normal distribution for kappa
%
% See: Gopal & Young,"von-Mises Fisher clusterig Models"
%

    D.penalizerparam = @penalizerparam;
    function penalizer_theta = penalizerparam(data)
        
        theta = estimatedefault(data);
        penalizer_theta.mu0 = theta.mu;
        penalizer_theta.kappa0 = theta.kappa;
        penalizer_theta.m = log(theta.kappa);
        penalizer_theta.sigma = 10;
        
    end

%% |sumparam|
% See <doc_distribution_common.html#18 distribution structure common members>.

    D.sumparam = @sumparam;
    function theta = sumparam(theta1, theta2)
        theta.mu = theta1.mu + theta2.mu;
        theta.kappa = theta1.kappa + theta2.kappa;
    end

%% |scaleparam|
% See <doc_distribution_common.html#19 distribution structure common members>.

    D.scaleparam = @scaleparam;
    function theta = scaleparam(scalar, theta)
        theta.mu = scalar * theta.mu;
        theta.kappa = scalar * theta.kappa;
    end

%% |sumgrad|
% See <doc_distribution_common.html#20 distribution structure common members>.

%% |scalegrad|
% See <doc_distribution_common.html#21 distribution structure common members>.

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
        str = [str, sprintf('kappa (scalar): %s\n', mat2str(theta.kappa, 4))];
        
        if nargout == 0
            str = [sprintf('%s distribution parameters:\n', D.name()), str];
        end
    end
    
%% |visualize|
%
% *Syntax*
%
%   handle_array = D.visualize(D, theta, vis_options)
%

    D.visualize = @visualize;
    function [varargout] = visualize(varargin)
        [varargout{1:nargout}] = vmf_visualize(D, varargin{:});
    end

%%

    D = mxe_addsharedfields(D);
end