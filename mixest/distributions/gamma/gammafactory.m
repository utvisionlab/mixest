%% |gammafactory|
% Construct a gamma distribution structure
%
% *Syntax*
%
%   D = gammafactory()
%
% *Description*
%
% |D = gammafactory()| returns a structure representing a gamma
% distribution.
%
% *Distribution Parameters*
%
% * *|a|* (positive scalar) : The shape parameter.
% * *|b|* (positive scalar) : The scale parameter.
%
% *Probability Density Function*
%
% The distribution has the following density:
% 
% $$ f(x)=\frac{1}{\Gamma(a)b^a} x^{a-1} \exp(-\frac{x}{b}) $$
%
% where $a > 0$ is the shape parameter, $b > 0$ is the scale parameter and
% $\Gamma$ is the gamma function.
% 
% *Example*
%
%   % Construct a gamma distribution
%   D = gammafactory();
%   % Plot the PDF for various values of the shape parameter:
%   data = 0:0.1:10;
%   f = zeros(4, numel(data));
%   for a = 1:4
%       theta = struct('a', a, 'b', 1);
%       f(a,:) = D.pdf(theta, data);
%   end
%   plot(data, f)
%   legend('a=1', 'a=2', 'a=3', 'a=4')
%
% <<img/gammafactory_01.png>>
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

function D = gammafactory(fixing)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'gamma';

%%

    thetamask = struct;
    fixed_a = false;
    fixed_b = false;
    
    % if the first argument is a struct with a 'fixing__' field, we should
    % make some parameters fixed
    if nargin>0 && isfield(fixing, 'fixing__')
        % fixing.thetamask is a parameter structure with non-existing or
        % NaN values for variable parameters and the desired fixed values
        % for fixed parameters.
        thetamask = fixing.thetamask;
        if isfield(thetamask, 'a') && ~isnan(thetamask.a)
            fixed_a = true;
        end
        if isfield(thetamask, 'b') && ~isnan(thetamask.b)
            fixed_b = true;
        end
    end

%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.

    D.M = build_param_manifold();
    function M = build_param_manifold()
        elements = struct;
        if ~fixed_a
            elements.a = positivefactory();
        end
        if ~fixed_b
            elements.b = positivefactory();
        end
        M = mxe_productmanifold(elements);
    end

%% |dim|
% See <doc_distribution_common.html#3 distribution structure common members>.

    D.dim = @dim; % parameter space dimensions
    function dim = dim()
        dim = 2;
        if fixed_a
            dim = dim - 1;
        end
        if fixed_b
            dim = dim - 1;
        end
    end

%% |datadim|
% See <doc_distribution_common.html#4 distribution structure common members>.

    D.datadim = @() 1; % data space dimensions

%%

    function store = ll_intermediate_params(data, store)
    % calculate intermediate parameters for ll 
        
        if ~isfield(store, 'sumLogData') || ~isfield(store, 'sumData') || ~isfield(store, 'sumW')
             data = mxe_readdata(data);
             weight = data.weight;
             N = data.size;
             data = data.data;
             
            if ~isfield(store, 'logData')
                store.LogData = log(data);
            end
            if ~isempty(weight)
                store.sumLogData = sum( store.LogData .* weight , 2 );
                store.sumData = sum ( data .* weight , 2 );
                store.sumW = sum(weight , 2);
            else
                store.sumLogData = sum( store.LogData , 2 );
                store.sumData = sum ( data , 2);    
                store.sumW = N;
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
       

        store = ll_intermediate_params(data, store);
        sumLogData = store.sumLogData;
        sumData = store.sumData;
        sumW = store.sumW;

        theta = fullparam(theta);
        ta = theta.a;
        tb = theta.b;
        
        ll = -sumW * (gammaln(ta) + ta*log(tb))  + (ta-1) * sumLogData ...
            - (1/tb) * sumData;
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
        
        theta = fullparam(theta);
        ta = theta.a;
        tb = theta.b;
        if ~isfield(store, 'logData')
            store.LogData = log(data);
        end
        llvec = -gammaln(ta) - ta*log(tb) + (ta-1)*store.LogData - (1/tb)*data;
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

        store = ll_intermediate_params(data, store);
        
        sumLogData = store.sumLogData;
        sumData = store.sumData;
        sumW = store.sumW;
        
        theta = fullparam(theta);
        ta = theta.a;
        tb = theta.b;
        
        if ~fixed_a
            % psi(0, ta) is digamma
            dll.a = -sumW *( log(tb) + psi(0, ta) ) + sumLogData;
        end
        if ~fixed_b
            dll.b = -sumW * ta/tb + (1/tb^2) * sumData;
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
        data = data.data;
        
        theta = fullparam(theta);
        ta = theta.a;
        tb = theta.b;
        
        dld = (ta-1)./data - 1/tb;
        if ~isempty(weight)
            dld = weight .* dld;
        end
    end

%% |cdf|
% See <doc_distribution_common.html#9 distribution structure common members>.

    D.cdf = @cdf;
    function y = cdf(theta, data)
        data = mxe_readdata(data);
        y = gammainc(data.data ./ theta.b, theta.a);
        %y = gamcdf(data.data, theta.a, theta.b); %matlab version
    end
    
%% |pdf|
% See <doc_distribution_common.html#10 distribution structure common members>.

%% |sample|
% See <doc_distribution_common.html#11 distribution structure common members>.

    if exist('randg','file') == 3
        % First check for matlab randg function
        samplegamma = @(x,y)randg(x, [1 y]);
    elseif exist('randgamma','file') == 3
        % Then Check for installed lightspeed randgamma
        samplegamma = @(x,y)randgamma(x * ones(1,y));
    else
        % At last run randraw which is twice slower
        samplegamma = @(x,y)randraw('gamma', x, [1 y]);
    end
    
    D.sample = @sample;
    function data = sample(theta, n)
        
        if nargin < 2, n = 1; end
        
        theta = fullparam(theta);
        data = samplegamma(theta.a, n);
        data = data * theta.b;

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
% Default estimation function for gamma distribution.
%
% *Syntax*
%
%   theta = D.estimatedefault(data)
%   theta = D.estimatedefault(data, options)
%   [theta, D] = D.estimatedefault(...)
%   [theta, D, info] = D.estimatedefault(...)
%   [theta, D, info, options] = D.estimatedefault(...)
%

%TODO doc

    D.estimatedefault = @estimatedefault;
    function [varargout] = estimatedefault(varargin)
        [varargout{1:nargout}] = gamma_estimatedefault(D, varargin{:}); 
    end

%% |penalizerparam|
% See <doc_distribution_common.html#15 distribution structure common members>.
%
% *Penalizer Info*
%
% The default penalizer is like having a exponential distribution prior on
% on parameter (a) and inversegamma on (a*b) as explained in the paper:
% Wiper et al. "Mixtures of Gamma Distributions with application"
%
% The exponential distribution has the following density:
% 
% $$ f(x) \propto \lambda \exp(-\lambda x) $$
%
% The inverse-gamma distribution has the following density:
%
% $$ f(x) \propto x^{-\alpha-1} \exp( -\beta / x ) $$
%

    D.penalizerparam = @penalizerparam;
    function penalizer_theta = penalizerparam(data)
        data = mxe_readdata(data);
        penalizer_theta.beta = mean(data.data)/2;
        penalizer_theta.alpha = 1.5;
        penalizer_theta.lambda = 0.01;
    end

%% |penalizercost|
% See <doc_distribution_common.html#16 distribution structure common members>.

    D.penalizercost = @penalizercost;
    function [costP, store] = penalizercost(theta, penalizer_theta, store)
        if nargin < 3
            store = struct;
        end
        theta = fullparam(theta);
        a = theta.a;
        b = theta.b;
        lambda = penalizer_theta.lambda;
        beta = penalizer_theta.beta;
        alpha = penalizer_theta.alpha;
        costP = - lambda * a;
        costP = costP - (alpha + 1) * ( log(a) + log(b) ) - beta / (b*a);
    end

%% |penalizergrad|
% See <doc_distribution_common.html#17 distribution structure common members>.

    D.penalizergrad = @penalizergrad;
    function [gradP, store] = penalizergrad(theta, penalizer_theta, store)
        if nargin < 3
            store = struct;
        end
        theta = fullparam(theta);
        a = theta.a;
        b = theta.b;
        lambda = penalizer_theta.lambda;
        beta = penalizer_theta.beta;
        alpha = penalizer_theta.alpha;
        if ~fixed_a
            gradP.a = -lambda - (alpha + 1)/a + beta / b / a.^2; 
        end
        if ~fixed_b
            gradP.b = - (alpha + 1)/b + beta / a / b.^2;
        end
    end

%% |sumparam|
% See <doc_distribution_common.html#18 distribution structure common members>.

    D.sumparam = @sumparam;
    function theta = sumparam(theta1, theta2)
        if ~fixed_a
            theta.a = theta1.a + theta2.a;
        end
        if ~fixed_b
            theta.b = theta1.b + theta2.b;
        end
    end

%% |scaleparam|
% See <doc_distribution_common.html#19 distribution structure common members>.

    D.scaleparam = @scaleparam;
    function theta = scaleparam(scalar, theta)
        if ~fixed_a
            theta.a = scalar * theta.a;
        end
        if ~fixed_b
            theta.b = scalar * theta.b;
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
        ta = theta.a;
        tb = theta.b;
        h = gammaln(ta) + ta + log(tb) + (1-ta) * psi(ta); 
    end

%% |kl|
% See <doc_distribution_common.html#23 distribution structure common members>.

    D.kl = @kl;
    function kl = kl(theta1, theta2)
        theta1 = fullparam(theta1);
        a1 = theta1.a;
        b1 = theta1.b;
        
        theta2 = fullparam(theta2);
        a2 = theta2.a;
        b2 = theta2.b;
        
        kl = (a1-a2).*psi(a1) + gammaln(a2) - gammaln(a1) + ...
            a2*(log(b2)-log(b1)) + a1/b2*b1 - a1;
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
        if fixed_a
            formatSpec = 'a: %g (fixed)\n';
        else
            formatSpec = 'a: %g\n';
        end
        str = [str, sprintf(formatSpec, theta.a)];
        
        if fixed_b
            formatSpec = 'b: %g (fixed)\n';
        else
            formatSpec = 'b: %g\n';
        end
        str = [str, sprintf(formatSpec, theta.b)];
        
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
        Afixing.thetamask = Athetamask;
        
        % construct the new distribution
        newD = gammafactory(Afixing);
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
        newD = gammafactory(Afixing);
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
        if fixed_a
            fulltheta.a = thetamask.a;
        end
        if fixed_b
            fulltheta.b = thetamask.b;
        end
    end

%%

    D = mxe_addsharedfields(D);
end