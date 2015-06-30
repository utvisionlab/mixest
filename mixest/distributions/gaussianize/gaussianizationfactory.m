%% |gaussianizationfactory|
% Construct a Gaussianization distribution structure
%
% *Syntax*
%
%   D = gaussianization(LayerD, num)
%   D = gaussianization(LayerD)
%
% *Description*
%
% |D = gaussianization(LayerD, num)| returns a structure representing a
% gaussianization distribution. |LayerD| is a distribution structure 
% defining gaussianization layer distribution type, and |num| is the number
% gaussianization layers
%
% |D = gaussianization(LayerD)| where |LayerD| is a cell array of
% distribution structures defined on the same data space, constructs a
% heterogeneous gaussianization distribution |D| where each layer 
% may be of a different distribution type.
%
% *Distribution Parameters*
%
% * *|D|* (|num-by-1| cell array of distribution parameter structures) :
% Contains the parameters for each layer.
%
% *Probability Density Function*
%
% The distribution is defined to be whitened Gaussian distribution in the
% transformed domain. The transformation is given by :
% 
% $$ y= G^{-1}(F_{num}(...G^{-1}(F_2(G^{-1}(F_1(x))))))) $$
%
% where $num$ is the number oflayers, $F_k$ represents the CDF of kth layer
% $G^{-1}$ is the inverse CDF of Gaussian.
% Note that $G^{-1}(F_k(.))$ is computer together as Gaussianization 
% for each layer
% 
% *Example*
%
%   % Construct a Gaussianization distribution of two layers
%   D = 
%   % Build a parameter structure for it:
%   theta{1} = 
%   theta{2} = 
%   % Plot the PDF:
%   x = 0:0.1:10;
%   y = 0:0.1:10;
%   surf
%
% <<img/...png>>
%
% For more detailes look at the following paper:
% Chen, S. S., & Gopinath, R. A. (2000). Gaussianization.

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Contributors:
%  Reshad Hosseini
%  Mohamadreza Mash'al
%
% Change log: 
%

function D = gaussianizationfactory(LayerD, num)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'gaussianization';

%%

    homogeneous = false;
    
    if nargin > 1
        homogeneous = true;
        Layers = num2cell(repmat(LayerD, [num,1]));
    else
        num = numel(LayerD);
        Layers = LayerD(:);
    end
     
%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.

    % Note: manifold and parameters are defined only on the variable layers
    if homogeneous
        D.M = powermanifold(LayerD.M, num);
    else
        elements = cell(num, 1);
        for it = 1:num
            elements{it} = Layers{it}.M;
        end
        D.M = mxe_productmanifold(elements);
    end



%% |num|
% Number of layers
%
% *Syntax*
%
%   num = D.num()
%
% *Description*
%
% |num = D.num()| returns the number of layers in the Gaussianization 
% distribution |D|. 
%
% *Note:* You need to include the parentheses |()| for this to work.
%

    D.num = @getnum;
    function result = getnum()
        result = num;
    end

%% |layer|
% layer distributions
%
% *Syntax*
%
%   D_k = D.layer(idx)
%
% *Description*
%
% |D_k = D.layer(idx)| where |idx| is the index of a layer, gives
% the layer distribution structure at |idx|.
%
% |D_k = D.layer(idx)| where |idx| is an index vector with more than
% one element, returns a cell array of the distribution structures indexed
% in |idx|.
%

    D.layer = @getlayer;
    function result = getlayer(idx)
        if numel(idx) == 1
            result = Layers{idx};
        else
            result = Layers(idx);
        end
    end


%% |addlayer|
% Add a Gaussianization layer (to the variable layers)
%
% *Syntax*
%
%   newD = D.addlayer(CmptD)
%   [newD, newtheta] = D.addlayer(CmptD, theta, CmptTheta)
%
% *Description*
%
% |newD = D.addlayer(CmptD)| returns |newD|, a Gaussianization distribution
% with the same layers as |D| plus the new variable layer
% distribution |CmptD|.
%
% |[newD, newtheta] = D.addlayer(CmptD, theta, CmptTheta)| also returns
% |newtheta|, the parameters for |newD|, given |theta|, the parameters for
% |D|, and |CmptTheta|, the parameters for |CmptD|.
%

    D.addlayer = @addlayer;
    function [newD, newtheta] = addlayer(CmptD, theta, CmptTheta)
        
        if nargout > 1
            newtheta = [theta(:); CmptTheta];
        end
        
        % build new Layers
        newLayers = [Layers; CmptD];
        
        % construct the new distribution
        newD = Gaussianizationfactory(newLayers);
    end

%% |removelayer|
% Remove a layer
%
% *Syntax*
%
%   newD = D.removelayer(idx)
%   [newD, newtheta] = D.removelayer(idx, theta)
%
% *Description*
%
% |newD = D.removelayer(idx)| returns |newD|, a Gaussianization distribution
% the same as |D| where the layer at index |idx| is removed from its
% layers.
%
% |[newD, newtheta] = D.removelayer(idx, theta)| also returns
% |newtheta|, the parameters for |newD|, given |theta|, the parameters for
% |D|.
%

    D.removelayer = @removelayer;
    function [newD, newtheta] = removelayer(idx, theta)
        
        if nargout > 1
            newtheta = theta;
            newtheta(idx) = [];
        end
        
        % build new Layers
        newLayers = Layers;
        newLayers(idx) = [];
        
        % construct the new distribution
        newD = Gaussianizationfactory(newLayers);
    end

%% |dim|
% See <doc_distribution_common.html#3 distribution structure common members>.

    D.dim = @dim; % parameter space dimensions
    function dim = dim()
        dim = 0; 
        for k = 1:num
            dim = dim + Layers{k}.dim(); % layer parameters
        end
    end

%% |datadim|
% See <doc_distribution_common.html#4 distribution structure common members>.

    D.datadim = @() Layers{1}.datadim(); % data space dimensions

%% |ll|
% See <doc_distribution_common.html#5 distribution structure common members>.

    D.ll = @ll;
    datadim = D.datadim();
    Gauss = mvnfactory(datadim);
    thetaG.mu = zeros(datadim,1);
    thetaG.sigma = eye(datadim);
    function [ll, store] = ll(theta, data, store)
        % layer by layer calculate difference between log-likelihood and
        % log-likelihood of Gaussian and sum over them

        store = struct;        
        data = mxe_readdata(data);
        ll_gauss = 0;
        ll = 0;
        for k = 1:num
            ll = ll + Layers{k}.ll(theta{k}, data) - ll_gauss;
            data = Layers{k}.gaussianize(theta{k}, data);
            ind = (data == inf);
            sind = sum(ind(:));
            if sind > 0
                warning('Some data become inf during Gaussianization');
                data(ind) = 3.17 + 0.01*randn(1,sind);
            end
            ind = (data == -inf);
            sind = sum(ind(:));
            if sind > 0
                warning('Some data become -inf during Gaussianization');
                data(ind) = -3.17 + 0.01*randn(1,sind);
            end
            % compute the likelihood of the data using Gaussian assumption
            ll_gauss = Gauss.ll(thetaG, data);
        end
        
    end

%% |llvec|
% See <doc_distribution_common.html#6 distribution structure common members>.

    D.llvec = @llvec;
    function [llvec, store] = llvec(theta, data, store)
        
        data = mxe_readdata(data);
        llvec_gauss = 0;
        llvec = 0;
        for k = 1:num
            llvec = llvec + Layers{k}.llvec(theta{k}, data) - llvec_gauss;
            data = Layers{k}.gaussianize(theta{k}, data);
            ind = (data == inf);
            sind = sum(ind(:));
            if sind > 0
                warning('Some data become inf during Gaussianization');
                data(ind) = 3.17 + 0.01*randn(1,sind);
            end
            ind = (data == -inf);
            sind = sum(ind(:));
            if sind > 0
                warning('Some data become -inf during Gaussianization');
                data(ind) = -3.17 + 0.01*randn(1,sind);
            end
            % compute the likelihood of the data using Gaussian assumption
            llvec_gauss = Gauss.llvec(thetaG, data);
        end
        
    end

%% |llgrad|
% See <doc_distribution_common.html#7 distribution structure common members>.

%% |llgraddata|
% See <doc_distribution_common.html#8 distribution structure common members>.

%% |cdf|
% See <doc_distribution_common.html#9 distribution structure common members>.

%% |pdf|
% See <doc_distribution_common.html#10 distribution structure common members>.

%% |sample|
% See <doc_distribution_common.html#11 distribution structure common members>.
%

%% |randparam|
% See <doc_distribution_common.html#12 distribution structure common members>.

%% |init|
% See <doc_distribution_common.html#13 distribution structure common members>.

    D.init = @init;
    function theta = init(data, varargin)

        data = mxe_readdata(data);
        data = data.data;
        theta = cell(num,1);
        for k = 1:num
            theta{k} = Layers{k}.init(data, varargin{:});
            data = Layers{k}.gaussianize(theta{k}, data);
        end
    end

%% |estimatedefault|
% Default estimation function for Gaussianization distribution. This function
% implements the expectation maximization (EM) method.
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
% This function supports the following options from the options described
% in <../estimation_options.html estimation options>.
%
% * *|theta0|*
% * *|verbosity|*
% * *|plotCost|*
% * *|crossVal|*
% * *|minIter|*
% * *|maxIter|*
% * *|maxTime|*
% * *|tolCost|*
% * *|tolCostDiff|*
% * *|statsfun|*
% * *|stopfun|*
%
% *Returned |info| fields*
%
% This function puts the following fields in the returned |info| structure
% array. You can read more about them in our documentation on
% <../estimation_statistics_structure.html estimation statistics
% structure>.
%
% * *|iter|*
% * *|cost|*
% * *|time|*
% * *|cvCost|*
%

    D.estimatedefault = @estimatedefault;
    function [varargout] = estimatedefault(data, varargin)
        if nargin > 1
            optionstot = mxe_options(varargin{1});
        else
            optionstot = mxe_options;
        end
        data = mxe_readdata(data);
        for k = 1:num
            options = optionstot;
            if ~isempty(optionstot.theta0)
                options.theta0 = optionstot.theta0{k};
            end
            if ~isempty(optionstot.penalizertheta)
                options.penalizertheta = optionstot.penalizertheta{k};
            end
            [varargout{1:nargout}] = Layers{k}.estimate(data, options, varargin{3:end});
            theta{k} = varargout{1};
            if nargout > 1
                Layersout{k} = varargout{2};
            end
            if nargout > 2
                infoout{k} = varargout{3};
            end
            if nargout > 3
                optionsout{k} = varargout{4};
            end
            data = Layers{k}.gaussianize(theta{k}, data);
            ind = (data == inf);
            sind = sum(ind(:));
            if sind > 0
                warning('Some data become inf during Gaussianization');
                data(ind) = 3.17 + 0.01*randn(1,sind);
            end
            ind = (data == -inf);
            sind = sum(ind(:));
            if sind > 0
                warning('Some data become -inf during Gaussianization');
                data(ind) = -3.17 + 0.01*randn(1,sind);
            end
        end
        varargout{1} = theta;
        if nargout > 1
            varargout{2} = gaussianizationfactory(Layersout);
        end
        if nargout > 2
            varargout{3} = infoout;
        end
        if nargout > 3
            varargout{4} = optionsout;
        end
    end

%% |penalizerparam|
% See <doc_distribution_common.html#15 distribution structure common members>.
%
% *Penalizer Info*
%
% The default penalizer for the Gaussianization distribution is the sum of the
% default penalizers of its layers.
%

    D.penalizerparam = @penalizerparam;
    function penalizer_theta = penalizerparam(data)
        penalizer_theta{1} = Layers{1}.penalizerparam(data);
        data = mxe_readdata(data);
        data = randn(size(data.data));
        for k = 1:num
            penalizer_theta{k} = Layers{k}.penalizerparam(data);
        end
    end

%% |penalizercost|
% See <doc_distribution_common.html#16 distribution structure common members>.

    D.penalizercost = @penalizercost;
    function [costP, store] = penalizercost(theta, penalizer_theta, store)
        
        if nargin < 3
            store = struct;
        end
        
        if ~isfield(store, 'layerStores')
            store.layerStores = cell(num,1);
        end
        
        costP = 0;
        for k = 1:num
            if ~isstruct(store.layerStores{k})
                store.layerStores{k} = struct;
            end
            [cost, store.layerStores{k}] = ...
                Layers{k}.penalizercost(theta{k}, penalizer_theta{k}, ...
                store.layerStores{k});
            costP = costP + cost;
        end
    end

%% |penalizergrad|
% See <doc_distribution_common.html#17 distribution structure common members>.


%% |sumparam|
% See <doc_distribution_common.html#18 distribution structure common members>.

    D.sumparam = @sumparam;
    function theta = sumparam(theta1, theta2)
        theta = cell(num, 1);
        for k = 1:num
            theta{k} = Layers{k}.sumparam(theta1{k}, theta2{k});
        end
    end

%% |scaleparam|
% See <doc_distribution_common.html#19 distribution structure common members>.

    D.scaleparam = @scaleparam;
    function theta = scaleparam(scalar, theta)
        for k = 1:num
            theta{k} = Layers{k}.scaleparam(scalar, theta{k});
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
        h = zeros(1,num);
        for k = 1:numtotal
            h(k) = Layers{k}.entropy(theta{k});
        end
    end

%% |kl|
% See <doc_distribution_common.html#23 distribution structure common members>.

    D.kl = @kl;
    function kl = kl(theta)
        kl = zeros(num);
        for k1 = 1:numtotal
            for k2 = 1:k1-1
                kl(k1,k2) = Layers{1}.kl(theta{k1}, theta{k2});
                kl(k2,k1) = kl(k1,k2);
            end
        end        
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
        for k = 1:num
            str = [str, sprintf('\nLayers(%d): %s \n', k, Layers{k}.name())]; %#ok<AGROW>
            str = [str, Layers{k}.display(theta{k})]; %#ok<AGROW>
        end
        
        if nargout == 0
            str = [sprintf('%s distribution parameters:\n', D.name()), str];
        end
    end

%% |visualize|
%


%%

    D = mxe_addsharedfields(D);
    
    % optimization options is for layer-wise optimization
    D.estimate = D.estimatedefault;
end