%% |mixturefactory|
% Construct a mixture distribution structure
%
% *Syntax*
%
%   D = mixturefactory(ComponentD, num)
%   D = mixturefactory(ComponentD)
%
% *Description*
%
% |D = mixturefactory(ComponentD, num)| returns a structure representing a
% mixture distribution. |ComponentD| is a distribution structure defining
% mixture component distribution type, and |num| is the number of mixture
% components.
%
% |D = mixturefactory(ComponentD)| where |ComponentD| is a cell array of
% distribution structures defined on the same data space, constructs a
% heterogeneous mixture distribution |D| where each component may be of a
% different distribution type.
%
% *Distribution Parameters*
%
% * *|D|* (|num-by-1| cell array of distribution parameter structures) :
% Contains the parameters for each component.
% * *|p|* (|num-by-1| vector) : The vector of component weights.
%
% *Probability Density Function*
%
% The distribution has the following density:
% 
% $$ f(x)=\sum_{k=1}^{num} \pi_k f_k(x) $$
%
% where $num$ is the number of components, $\pi_k$ is the weight for the
% k-th component, and $f_k$ represents the density function for the k-th
% component.
% 
% *Example 1*
%
%   % Construct a mixture of three bivariate normal distributions:
%   D = mixturefactory(mvnfactory(2), 3);
%   % Build a parameter structure for it:
%   theta.D{1} = struct('mu', [3; -1], 'sigma', [2 1; 3 4]);
%   theta.D{2} = struct('mu', [0; 3], 'sigma', [3 0; 0 4]);
%   theta.D{3} = struct('mu', [-4; -1], 'sigma', [4 0; 0 3]);
%   theta.p = [0.4; 0.3; 0.3];
%   % Plot the PDF:
%   x = -10:0.2:10;
%   y = -10:0.2:10;
%   [X, Y] = meshgrid(x, y);
%   data = [X(:) Y(:)]';
%   f = D.pdf(theta, data);
%   surf(X, Y, reshape(f, size(X)));
% 
% <<img/mixturefactory_01.png>>
% 
% *Example 2*
%
%   % Construct a heterogeneous mixture distribution
%   % containing a gamma and a Gaussian component:
%   D = mixturefactory({gammafactory(); mvnfactory(1)});
%   % Build a parameter structure for it:
%   theta.D{1} = struct('a', 2, 'b', 1);
%   theta.D{2} = struct('mu', 7, 'sigma', 4);
%   theta.p = [0.7; 0.3];
%   % Plot the PDF:
%   x = 0:0.1:10;
%   plot(x, D.pdf(theta, x))
%
% <<img/mixturefactory_02.png>>
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

function D = mixturefactory(ComponentD, num)

%% |name|
% See <doc_distribution_common.html#1 distribution structure common members>.

    D.name = @() 'mixture';
    
%% 
% Flag to control the memory usage (resulting code will be slower)

    low_memory = false;
    
%%

    homogeneous = false;
    hXcache = []; % ll cache for fixed components (on train data)
    cacheValid = false; % flag indicating whether ll-cache is up-to-date
    
    % if the first argument is a struct with a 'fixing__' field, we should
    % make some parameters fixed
    fixing = ComponentD;
    if nargin>0 && isfield(fixing, 'fixing__')
        
        varfixupdate(fixing.varD, fixing.fixedD, fixing.fixedtheta, true);
        
        if isfield(fixing, 'hXcache')
            hXcache = fixing.hXcache;
            cacheValid = true;
        end
        if isfield(fixing, 'data') && isfield(fixing, 'datapatchsize')
            fill_cache(fixing.data, fixing.datapatchsize);
        elseif isfield(fixing, 'data')
            fill_cache(fixing.data);
        end
        
    else



        if nargin > 1
            homogeneous = true;
            Components = num2cell(repmat(ComponentD, [num,1]));
        else
            num = numel(ComponentD);
            Components = ComponentD(:);
        end
        % variables used for partial estimation (these are managed in the varfixupdate function)
        varD = Components;
        fixedD = {};
        fixedtheta = struct('D', {}, 'p', {});
        numfixed = 0; % number of fixed components. Determinig whether we have fixed components (and that theta.p has an extra element) should be done by checking this.
        numtotal = num; % numtotal = num + numfixed; is set every time num/numfixed is changed (because it is used a lot)
        nump = num; % number of elements in theta.p (you may use this instead of checking numfixed)
    end

%%

    function varfixupdate(AvarD, AfixedD, Afixedtheta, invalidateCache)
    % performs updates required after changing var/fixed components
    % Note: After calling this, you should update the cache by calling fillcache(data), if you have data
        
        if nargin < 4
            invalidateCache = true;
        end
        
        varD = AvarD(:);
        fixedD = AfixedD(:);
        num = numel(varD);
        numfixed = numel(fixedD);
        numtotal = num + numfixed;
        Components = [varD; fixedD];
        if numfixed > 0
            nump = num + 1;
        else
            nump = num;
        end
        
        % fix up Afixedtheta
        if numfixed > 0
            if isempty(Afixedtheta)
                % empty fixedTheta => generate a random fixedTheta
                Afixedtheta = struct;
                Afixedtheta.D = cell(numfixed, 1);
                for k = 1:numfixed
                    Afixedtheta.D{k} = fixedD{k}.randparam();
                end
                Afixedtheta.p = ones(numfixed, 1) ./ numfixed;
            elseif iscell(Afixedtheta)
                % fixedTheta is a cell array containing parameters for each fixed component
                temp = Afixedtheta(:);
                Afixedtheta = struct;
                Afixedtheta.D = temp; clear temp        
                Afixedtheta.p = ones(numfixed, 1) ./ numfixed;
            else
                if isfield(Afixedtheta, 'D')
                    % fixedTheta is truly a mixture parameter
                    Afixedtheta.D = Afixedtheta.D(:);
                    Afixedtheta.p = Afixedtheta.p(:);                
                    if numel(Afixedtheta.p) == nump
                        % if e.g. fixedTheta comes from a call to parmix.subparam()
                        % and contains an additional weight, remove it.
                        Afixedtheta.p = Afixedtheta.p(1:numfixed);
                    end
                    % fixedTheta.p must always sum to 1
                    Afixedtheta.p = Afixedtheta.p ./ sum(Afixedtheta.p);
                else
                    % parameters for a single component is given
                    Afixedtheta.D = {Afixedtheta};
                    Afixedtheta.p = 1;                
                end
            end
        else
            Afixedtheta = struct('D', {}, 'p', {});
        end
        fixedtheta = Afixedtheta;
        
%         % update the parameter manifold
%         D.M = build_param_manifold();
%         % We need to update functions added by mxe_addsharedfields, when D.M is changed
%         D = mxe_addsharedfields(D);
        
        % invalidate the ll cache
        if invalidateCache
            cacheValid = false;
        end
    end

%%

    function fill_cache(all_data, datapatchsize)
    % caches ll for fixed components on all data (including
    % cross-validation data, other patches, etc).
    %
    % Note: (numfixed, fixedD, fixedtheta) must be set before calling this
        
        if nargin < 2
            datapatchsize = inf;
        end
        hXcache = calc_cache(all_data, datapatchsize);
        cacheValid = true;
    end

    function hX = calc_cache(data, datapatchsize)
    % calculate ll-cache for given data (used if cache is not filled)
    %
    % Note: (numfixed, fixedD, fixedtheta) must be set before calling this
        
        % copied from weighting
        data = mxe_readdata(data, false);
        
        hX = zeros(numfixed, data.size);
        if numfixed == 0
            return
        end
        
        data.weight = [];
        %data.index = []; % calculate the cache on all data  %TODO patching
        
        if isinf(datapatchsize)
            datapatchsize = data.size;
        end
        
        numk2 = ceil(data.size / datapatchsize);
        
        for k2 = 1:numk2
            
            bind = 1 + (k2-1)*datapatchsize;
            eind = min(k2*datapatchsize, data.size);
            data.index = bind:eind; % done with patching
            for k = 1:numfixed
                hX(k,bind:eind) = log(fixedtheta.p(k)) + ...
                    fixedD{k}.llvec(fixedtheta.D{k}, data);
            end 
            
        end
        % to minimize memory use, we sum up the fixed lls into one row.
        % (note: we can safely remove this line)
        hX = logsumexp(hX, 1);
    end

%% |M|
% See <doc_distribution_common.html#2 distribution structure common members>.

    % Note: manifold and parameters are defined only on the variable components
    D.M = build_param_manifold();
    function M = build_param_manifold()
        
        if num == 0
            M = struct;
            return
        end
        
        if homogeneous
            varM = powermanifold(ComponentD.M, num);
        else
            elements = cell(num, 1);
            for k = 1:num
                elements{k} = varD{k}.M;
            end
            varM = mxe_productmanifold(elements);
        end
        % When having fixed components, p(num+1) is the global weight for fixedtheta.p
        % Note: fixedtheta.p always sums to 1
        p = simplexfactory(nump);
        M = productmanifold(struct('D', varM, 'p', p));
    end

    % from here on use Components{k} to access k'th component irrespective
    % to the mixture being homogeneous or heterogeneous.


%% |num|
% Number of components (excluding any fixed components)
%
% *Syntax*
%
%   num = D.num()
%
% *Description*
%
% |num = D.num()| returns the number of components in the mixture |D|. When
% some components are fixed, only the number of variable components is
% returned.
%
% *Note:* You need to include the parentheses |()| for this to work.
%

    D.num = @getnum;
    function result = getnum()
        result = num;
    end
    
%% |numfixed|
% Number of fixed components
%
% *Syntax*
%
%   numfixed = D.numfixed()
%
% *Description*
%
% |numfixed = D.numfixed()| returns the number of fixed components in the
% mixture |D|.
%
% *Note:* You need to include the parentheses |()| for this to work.
%

    D.numfixed = @getnumfixed;
    function result = getnumfixed()
        result = numfixed;
    end

%% |numtotal|
% Total number of components
%
% *Syntax*
%
%   numtotal = D.numtotal()
%
% *Description*
%
% |numtotal = D.numtotal()| returns the total number of components
% (variable and fixed) in the mixture |D|.
%
% *Note:* You need to include the parentheses |()| for this to work.
%

    D.numtotal = @getnumtotal;
    function result = getnumtotal()
        result = numtotal;
    end

%% |nump|
% Number of elements in the component weight vector (theta.p)
%
% *Syntax*
%
%   nump = D.nump()
%
% *Description*
%
% |nump = D.nump()| returns the number of elements in the component weight
% vector in the parameters for the mixture |D|. When there are no fixed
% components, this is the same as the number of mixture components. When
% some components are fixed, |nump| equals the number of variable
% components plus one, since in addition to the weights for the variable
% components, we store another weight which can scale the weights of the
% fixed components uniformly.
%
% *Note:* You need to include the parentheses |()| for this to work.
%

    D.nump = @getnump;
    function result = getnump()
        result = nump;
    end

%% |component|
% Component distributions
%
% *Syntax*
%
%   D_k = D.component(idx)
%
% *Description*
%
% |D_k = D.component(idx)| where |idx| is the index of a component, gives
% the component distribution structure at |idx|.
%
% |D_k = D.component(idx)| where |idx| is an index vector with more than
% one element, returns a cell array of the distribution structures indexed
% in |idx|.
%
% Valid range for the indices is from 1 to |D.numtotal()|.
%

    D.component = @getcomponent;
    function result = getcomponent(idx)
        if numel(idx) == 1
            result = Components{idx};
        else
            result = Components(idx);
        end
    end

%% |varD|
% Variable component distributions
%
% *Syntax*
%
%   varD_k = D.varD(idx)
%
% *Description*
%
% |varD_k = D.varD(idx)| where |idx| is the index of a variable component,
% gives the variable component distribution structure at |idx|.
%
% |varD_k = D.varD(idx)| where |idx| is an index vector with more than one
% element, returns a cell array of the variable distribution structures
% indexed in |idx|.
%
% Valid range for the indices is from 1 to |D.num()|.
%

    D.varD = @getvarD;
    function result = getvarD(idx)
        if numel(idx) == 1
            result = varD{idx};
        else
            result = varD(idx);
        end
    end

%% |fixedD|
% Fixed component distributions
%
% *Syntax*
%
%   fixedD_k = D.fixedD(idx)
%
% *Description*
%
% |fixedD_k = D.fixedD(idx)| where |idx| is the index of a fixed component,
% gives the fixed component distribution structure at |idx|.
%
% |fixedD_k = D.fixedD(idx)| where |idx| is an index vector with more than
% one element, returns a cell array of the fixed distribution structures
% indexed in |idx|.
%
% Valid range for the indices is from 1 to |D.numfixed()|.
%

    D.fixedD = @getfixedD;
    function result = getfixedD(idx)
        if numel(idx) == 1
            result = fixedD{idx};
        else
            result = fixedD(idx);
        end
    end

%% |fixedparam|
% Parameters stored for the fixed components
%
% *Syntax*
%
%   fixedtheta = D.fixedparam()
%
% *Description*
%
% |fixedtheta = D.fixedparam()| returns the parameter structure stored for
% the fixed components in the mixture |D|. The number of elements in
% |fixedtheta.D| and |fixedtheta.p| equals |D.numfixed()|.
%
% *Note:* You need to include the parentheses |()| for this to work.
%

    D.fixedparam = @getfixedparam;
    function result = getfixedparam()
        result = fixedtheta;
    end

%% |subparam|
% Extract a subset of component parameters
%
% *Syntax*
%
%   subtheta = D.subparam(theta, idx)
%
% *Description*
%
% |subtheta = D.subparam(theta, idx)| returns the mixture parameter
% structure |subtheta| containing the parameters for the subset of
% components indexed by |idx|, from |theta|. An additional weight is added
% to |subtheta.p| to make the sum of the weights equal to one.
%
% When some components are fixed, you can pass indices from |D.num() + 1|
% to |D.num() + D.numfixed()|, to refer to fixed components.
%

    D.subparam = @subparam;
    function subtheta = subparam(theta, idx)
    % get a subset of parameters from theta, that are indicated by
    % component indices in idx.
    
        if any(idx > num) && numfixed>0
            theta = fullparam(theta);
        end
        subtheta.D = theta.D(idx);
        theta.p = theta.p(:);
        newp = theta.p(idx);
        w = 1 - sum(newp);
        subtheta.p = [newp; w];
    end

%% |fullparam|
% Get the parameters for the entire mixture, given variable component
% parameters
%
% *Syntax*
%
%   fulltheta = D.fullparam(theta)
%
% *Description*
%
% |fulltheta = D.fullparam(theta)| can be used when some components are
% fixed, to obtain the parameters for all the components, given the
% partial parameters |theta| corresponding to the variable components.
%

    D.fullparam = @fullparam;
    function fulltheta = fullparam(theta)
    % returns parameters for the entire mixture distribution (including
    % variable and fixed component parameters)
    
        if numel(theta.D) == numtotal
            fulltheta = theta;
            return
        end
        fulltheta.D = [theta.D(:); fixedtheta.D];
        w = theta.p(num+1);
        theta.p = theta.p(1:num);
        fulltheta.p = [theta.p(:); w.*fixedtheta.p];
    end

%% |addcomponent|
% Add a mixture component (to the variable components)
%
% *Syntax*
%
%   newD = D.addcomponent(CmptD)
%   [newD, newtheta] = D.addcomponent(CmptD, theta, CmptTheta)
%
% *Description*
%
% |newD = D.addcomponent(CmptD)| returns |newD|, a mixture distribution
% with the same components as |D| plus the new variable component
% distribution |CmptD|.
%
% |[newD, newtheta] = D.addcomponent(CmptD, theta, CmptTheta)| also returns
% |newtheta|, the parameters for |newD|, given |theta|, the parameters for
% |D|, and |CmptTheta|, the parameters for |CmptD|.
%

    D.addcomponent = @addcomponent;
    function [newD, newtheta] = addcomponent(CmptD, theta, CmptTheta)
        
        if nargout > 1
            newtheta.D = [theta.D(:); CmptTheta.D];
            newtheta.p = [theta.p(:); CmptTheta.p];
        end
        
        % build new varD
        newVarD = [varD; CmptD];
        
        % build fixing struct
        Afixing = struct('fixing__', true);
        Afixing.varD = newVarD;
        Afixing.fixedD = fixedD;
        Afixing.fixedtheta = fixedtheta;
        if cacheValid
            Afixing.hXcache = hXcache; % cache is still valid
        end
        
        % construct the new distribution
        newD = mixturefactory(Afixing);
    end

%% |removecomponent|
% Remove a component
%
% *Syntax*
%
%   newD = D.removecomponent(idx)
%   [newD, newtheta] = D.removecomponent(idx, theta)
%
% *Description*
%
% |newD = D.removecomponent(idx)| returns |newD|, a mixture distribution
% the same as |D| where the component at index |idx| is removed from its
% components.
%
% |[newD, newtheta] = D.removecomponent(idx, theta)| also returns
% |newtheta|, the parameters for |newD|, given |theta|, the parameters for
% |D|.
%

    D.removecomponent = @removecomponent;
    function [newD, newtheta] = removecomponent(idx, theta)
        
        if nargout > 1
            newtheta = theta;
            newtheta.D(idx) = [];
            newtheta.p(idx) = [];
        end
        
        % build new varD
        newVarD = varD;
        newVarD(idx) = [];
        
        % build fixing struct
        Afixing = struct('fixing__', true);
        Afixing.varD = newVarD;
        Afixing.fixedD = fixedD;
        Afixing.fixedtheta = fixedtheta;
        if cacheValid
            Afixing.hXcache = hXcache; % cache is still valid
        end
        
        % construct the new distribution
        newD = mixturefactory(Afixing);
    end

%% |invertindex|
% Obtain an inverted component index vector (Find the other components)
%
% *Syntax*
%
%   invidx = D.invertindex(idx)
%   invidx = D.invertindex(idx, 'fixed')
%
% *Description*
%
% |invidx = D.invertindex(idx)| returns the indices of the variable
% components other than those indexed  in |idx|.
%
% |invidx = D.invertindex(idx, 'fixed')| returns the indices of the fixed
% components other than those indexed  in |idx|.
%
% *Example*
%
%   D = mixturefactory(mvnfactory(1), 6);
%   invidx = D.invertindex([3,5])
%   
%   invidx =
%      1     2     4     6
%

    D.invertindex = @invertindex;
    function invidx = invertindex(idx, fixedflag)
        
        if nargin<2 || strcmpi(fixedflag, 'var')
            
            invidx = 1:num;
            invidx(idx) = [];
            
        elseif strcmpi(fixedflag, 'fixed')
            
            invidx = 1:numfixed;
            invidx(idx) = [];
            
        else
            error('mixture.invertindex: syntax error')
        end
    end

%% |fixate|
% Make some variable component(s) fixed
%
% *Syntax*
%
%   newD = D.fixate(idx, theta)
%   [newD, theta] = D.fixate(idx, theta)
%   [newD, theta] = D.fixate(idx, theta, data)
%   [newD, theta, idxFixed] = D.fixate(...)
%   [newD, theta, idxFixed, idxMap] = D.fixate(...)
%
% *Description*
%
% |newD = D.fixate(idx, theta)| returns |newD|, the same mixture as |D|
% where the variable components indexed by the index vector |idx| are fixed
% and their parameter values will not change during estimation. |theta| is
% a full parameter structure for the mixture. The values for the parameters
% of the newly fixed components are read from |theta|.
%
% |[newD, theta] = D.fixate(idx, theta)| also returns the output |theta|,
% containing the parameters for the variable components extracted from the
% input |theta|, applicable to the new mixture, |newD|.
%
% |[newD, theta] = D.fixate(idx, theta, data)| can be used to fill the
% log-likelihood cache for the fixed components, using the given |data|.
% The cache is used for better performance while calculating the
% log-likelihood for |newD|.
%
% |[newD, theta, idxFixed] = D.fixate(...)| also returns |idxFixed|, the
% indices of the newly fixed components in |newD.fixedD()|.
%
% |[newD, theta, idxFixed, idxMap] = D.fixate(...)| also returns |idxMap|,
% an index-map vector, with length |D.num()|, which maps the old indices of
% the variable components in |D| to their new indices in |newD| (for those
% that are not being fixed). You can use the map like this: |idx_new =
% idxMap(idx_old)|, where |idx_old| denotes an index for a variable
% component in |D| and |idx_new| denotes its new index in the variable
% components of |newD|.
%
% *Example*
%
% This is how the |estimatepartial| function uses |invertindex|, |fixate|,
% |estimate| and |unfix| to perform partial parameter estimation for
% components indicated by |idx|:
%
%   function [newtheta, newD, info, options] = estimatepartial(idx, theta, data, options)
%     % make the other components fixed
%     invidx = invertindex(idx);
%     [newD, theta0, idxfixed] = fixate(invidx, theta, data);
%     % use the given theta as the initial point for estimation
%     options.theta0 = theta0;
%     % estimate the resulting partial mixture parameters
%     [newtheta, newD, info, options] = newD.estimate(data, options);
%     % bring the fixed components back
%     [newD, newtheta] = newD.unfix(idxfixed, invidx, newtheta);
%   end
%

    D.fixate = @fixate;
    function [newD, theta, idxFixed, idxMap] = fixate(idx, theta, data, datapatchsize)

        % validation
        if num == 0
            newD = D;
            idxFixed = [];
            idxMap = 1:num;
            return
        end
        
        if strcmpi(idx, 'all')
            idx = 1:num;
        else
            if any(idx > num)
                idx(idx > num) = [];
                warning('mixture.fixate: Indices more than num were removed')
            end
        end
        if isempty(idx)
            newD = D;
            idxFixed = [];
            idxMap = 1:num;
            return
        end
        
        theta.D = theta.D(:);
        theta.p = theta.p(:);

        
        % compute the new mixture structure
        if numfixed > 0
            newFixedD = [fixedD; varD(idx)];
            newFixedTheta.D = [fixedtheta.D; theta.D(idx)];
            w = theta.p(num+1);
            newFixedTheta.p = [w.*fixedtheta.p; theta.p(idx)]; % will be normalized
        else
            newFixedD = varD(idx);
            newFixedTheta.D = theta.D(idx);
            newFixedTheta.p = theta.p(idx); % will be normalized
        end
        newVarD = varD;
        newVarD(idx) = [];
        
        
        % output new indices
        if nargout > 2
            idxFixed = 1:numel(idx);
            idxFixed = idxFixed + numfixed;
        end
        
        % calculate the index map for original (variable) components.
        % can be used like:   idx_new = idxMap(idx_old)
        if nargout > 3
            idxMap = zeros(1,num);
            invidx = invertindex(idx);
            if numel(invidx) > 0
                idxMap(invidx) = 1:numel(invidx);
            end
        end
        
        % update the variable theta
        if nargout > 1
            theta.D(idx) = [];
            newp = theta.p(1:num);
            newp(idx) = [];
            w = 1 - sum(newp);
            theta.p = [newp; w];
        end
        
        % build fixing struct
        Afixing = struct('fixing__', true);
        Afixing.varD = newVarD;
        Afixing.fixedD = newFixedD;
        Afixing.fixedtheta = newFixedTheta;
        if nargin > 2
            Afixing.data = data;
        end
        if nargin > 3
            Afixing.datapatchsize = datapatchsize;
        end
        
        % construct the new distribution
        newD = mixturefactory(Afixing);
    end

%% |unfix|
% Make some fixed component(s) variable (undo |fixate|)
%
% *Syntax*
%
%   newD = D.unfix(idx)
%   [newD, theta] = D.unfix(idx, idxFixated, theta)
%   [newD, theta] = D.unfix(idx, idxFixated, theta, data)
%   [newD, theta, idxUnfixed] = D.unfix(...)
%   [newD, theta, idxUnfixed, idxMap] = D.unfix(...)
%
% *Description*
%
% |newD = D.unfix(idx)| returns |newD|, the same mixture as |D| where the
% fixed components indexed by the index vector |idx| are made variable for
% estimation.
%
% |[newD, theta] = D.unfix(idx, idxFixated, theta)| brings back the fixed
% components to their original place before fixation. |idxFixated|, same
% size as |idx|, should be the index vector that was previously given to
% |fixate| (as its first argument) to fixate these components. When
% |idxFixated| is an empty vector, the components are added to the end of
% the variable components. The input |theta| is a parameter structure for
% the mixture |D| which has some fixed components, and the output |theta|
% is the parameter structure for the mixture |newD| where the indexed
% components are made variable.
%
% |[newD, theta] = D.unfix(idx, idxFixated, theta, data)| can be used to
% fill the log-likelihood cache for any remaining fixed components, using
% the given |data|. The cache is used for better performance while
% calculating the log-likelihood for |newD|.
%
% |[newD, theta, idxUnfixed] = D.unfix(...)| also returns |idxUnfixed|, the
% indices of the unfixed components in |newD.varD()|.
%
% |[newD, theta, idxUnfixed, idxMap] = D.unfix(...)| also returns |idxMap|,
% an index-map vector, with length |D.numfixed()|, which maps the old
% indices of the fixed components in |D| to their new indices in |newD|
% (for those that are not being made variable). You can use the map like
% this: |idx_new = idxMap(idx_old)|, where |idx_old| denotes an index for a
% fixed component in |D| and |idx_new| denotes its new index in the fixed
% components of |newD|.
%
% *Example*
%
% This is how the |estimatepartial| function uses |invertindex|, |fixate|,
% |estimate| and |unfix| to perform partial parameter estimation for
% components indicated by |idx|:
%
%   function [newtheta, newD, info, options] = estimatepartial(idx, theta, data, options)
%     % make the other components fixed
%     invidx = invertindex(idx);
%     [newD, theta0, idxfixed] = fixate(invidx, theta, data);
%     % use the given theta as the initial point for estimation
%     options.theta0 = theta0;
%     % estimate the resulting partial mixture parameters
%     [newtheta, newD, info, options] = newD.estimate(data, options);
%     % bring the fixed components back
%     [newD, newtheta] = newD.unfix(idxfixed, invidx, newtheta);
%   end
%

    D.unfix = @unfix;
    function [newD, theta, idxUnfixed, idxMap] = unfix(idx, idxFixated, theta, data)
    % idxFixated: idx in variable components (that was previously given to
    % the fixate function) to be reverted, or [] to append to the last of
    % variable components
        
        % validation
        if numfixed == 0
            newD = D;
            idxUnfixed = [];
            idxMap = 1:numfixed;
            return
        end
        
        if strcmpi(idx, 'all')
            idx = 1:numfixed;
        else
            if any(idx > numfixed)
                idx(idx > numfixed) = [];
                warning('mixture.unfix: Indices more than numfixed were removed')
            end
        end
        if isempty(idx)
            newD = D;
            idxUnfixed = [];
            idxMap = 1:numfixed;
            return
        end

        
        % calculate required indices in new variable components
        if nargin<1 || isempty(idxFixated)
            % no idxFixated given => append to the end
            idxFixated = 1:numel(idx);
            idxFixated = idxFixated + num;
        end
        newNum = num + numel(idx); % new number of variable components
        invIdxFixated = 1:newNum;
        invIdxFixated(idxFixated) = [];
        
        % update the variable theta
        if nargout > 1
            tempD = theta.D(:);
            theta.D = cell(newNum, 1);
            theta.D(idxFixated) = fixedtheta.D(idx);
            theta.D(invIdxFixated) = tempD;
            tempp = theta.p(:);
            w = tempp(num+1);
            theta.p = zeros(newNum, 1);
            theta.p(idxFixated) = w.*fixedtheta.p(idx);
            theta.p(invIdxFixated) = tempp(1:num);
            % if there will remain other fixed components, add the extra weight to theta.p
            if numfixed-numel(idx) > 0
                w = 1 - sum(theta.p);
                theta.p = [theta.p; w];
            end
        end
        
        % compute the new mixture structure
        newVarD = cell(newNum, 1);
        newVarD(idxFixated) = fixedD(idx);
        newVarD(invIdxFixated) = varD;
        newFixedD = fixedD;
        newFixedD(idx) = [];
        newFixedTheta = fixedtheta;
        newFixedTheta.D(idx) = [];
        newFixedTheta.p(idx) = []; % will be normalized
        
        
        % output new indices
        if nargout > 2
            idxUnfixed = idxFixated;
        end
        
        % calculate the index map for original fixed components. 
        % can be used like:   idx_new = idxMap(idx_old)
        if nargout > 3
            idxMap = zeros(1,numfixed);
            invidx = invertindex(idx, 'fixed');
            if numel(invidx) > 0
                idxMap(invidx) = 1:numel(invidx);
            end
        end
        
        % build fixing struct
        Afixing = struct('fixing__', true);
        Afixing.varD = newVarD;
        Afixing.fixedD = newFixedD;
        Afixing.fixedtheta = newFixedTheta;
        if nargin > 3
            Afixing.data = data;
        end
        
        % construct the new distribution
        newD = mixturefactory(Afixing);
    end

%% |splitinit|
% Calculate the initial parameters for two splitted mixture components to
% be substituted for the given component.
%
% *Syntax*
%
%   newtheta = D.splitinit(idx, theta)
%   newtheta = D.splitinit(idx, theta, options)
%   newtheta = D.splitinit(idx, theta, options, data)
%   [newtheta, store] = D.splitinit(idx, theta, options, data, store)
%
% *Description*
%
% |newtheta = D.splitinit(idx, theta)| returns the parameters for a new
% mixture resulted by splitting the component indexed by |idx|. The input
% |theta| is the parameters for the mixture before the split.
% Initialization for the splitted components is performed using the default
% methods for each parameter.
%
% |newtheta = D.splitinit(idx, theta, options)| uses |options.splitInit| to
% initialize the parameters for the splitted components as described
% <../estimation_options.html#10 here>.
%
% |newtheta = D.splitinit(idx, theta, options, data)| can be used if the
% initialization method requires estimation data.
%
% |[newtheta, store] = D.splitinit(idx, theta, options, data, store)| can
% be used for caching purposes as described <../caching.html here>.
%
% *Note:* This function is automatically called by the |split| function and
% you don't usually need to call this function separately.
%

    D.splitinit = @splitinit;
    function [varargout] = splitinit(varargin)
        [varargout{1:nargout}] = mixture_splitinit(D, varargin{:}); 
    end

%% |split|
% Split a variable component in two
%
% *Syntax*
%
%   newD = D.split(idx)
%   [newD, newtheta] = D.split(idx, theta)
%   [newD, newtheta] = D.split(idx, theta, options)
%   [newD, newtheta] = D.split(idx, theta, options, data)
%   [newD, newtheta, idxSplitted, idxMap] = D.split(...)
%   [newD, newtheta, idxSplitted, idxMap, store] = D.split(idx, theta, options, data, store)
%
% *Description*
%
% |newD = D.split(idx)| returns the new mixture, |newD|, resulted by
% splitting the component indexed by |idx| from |D|.
%
% |[newD, newtheta] = D.split(idx, theta)| also initializes the parameters
% for the splitted components. The input |theta| is the parameter structure
% for the mixture |D| before the split. The output |newtheta| is the
% parameter structure of the new mixture, |newD|, after the split
% (containing the parameters for an additional component). The initialized
% parameter values for the first splitted component are stored at
% |newtheta.D{idx}| and for the second at |newtheta.D{end}|. The component
% weights in |newtheta.p| are also updated accordingly, such that each
% splitted component has half the weight of the original component.
%
% |[newD, newtheta] = D.split(idx, theta, options)| uses
% |options.splitInit| to initialize the parameters for the splitted
% components as described <../estimation_options.html#10 here>.
%
% |[newD, newtheta] = D.split(idx, theta, options, data)| can be used if
% the initialization method requires estimation data.
%
% |[newD, newtheta, idxSplitted, idxMap] = D.split(...)| also returns the
% two indices of the splitted components in the vector |idxSplitted|. The
% output |idxMap| is an index-map vector, with length |D.num()|, which maps
% the old indices of the components in |D| to their new indices in |newD|.
% You can use the map like this: |idx_new = idxMap(idx_old)|, where
% |idx_old| denotes an index for a component in |D| and |idx_new| denotes
% its new index in the components of |newD|.
%
% |[newD, newtheta, idxSplitted, idxMap, store] = D.split(idx, theta,
% options, data, store)| can be used for caching purposes as described
% <../caching.html here>.
%
% *Example*
%
% The following code demonstrates a basic usage of split-related functions:
%
%   % Construct a mixture distribution and generate random 
%   % parameters and data for demonstration
%   D = mixturefactory(mvnfactory(2), 3);
%   theta = D.randparam();
%   data = D.sample(theta, 1000);
%   % Find candidate components for splitting
%   idx_split = D.splitcandidates(theta, data);
%   % Split the first candidate component
%   idx = idx_split(1);
%   [newD, newtheta, idxSplitted] = D.split(idx, theta);
%   % We may estimate the parameters merely for the splitted components
%   newtheta = newD.estimatepartial(idxSplitted, newtheta, data)
%

    D.split = @split;
    function [newD, newtheta, idxSplitted, idxMap, store] = split(idx, theta, options, data, store)
        
        if nargout > 1
            if nargin < 5
                store = struct;
            end
            if nargin < 4
                data = [];
            end
            if nargin < 3
                options = mxe_options();
            else
                options = mxe_options(options);
            end
            [newtheta, store] = splitinit(idx, theta, options, data, store);
        end
        
        % calculate new indices
        if nargout > 2
            idxSplitted = [idx, num+1];
        end
        
        % calculate the index map for original (variable) components.
        % can be used like:   idx_new = idxMap(idx_old)
        if nargout > 3
            idxMap = 1:num;
        end
        
        % build new varD
        newVarD = [varD; varD(idx)];
        
        % build fixing struct
        Afixing = struct('fixing__', true);
        Afixing.varD = newVarD;
        Afixing.fixedD = fixedD;
        Afixing.fixedtheta = fixedtheta;
        if cacheValid
            Afixing.hXcache = hXcache; % cache is still valid
        end
        
        % construct the new distribution
        newD = mixturefactory(Afixing);
    end

%% |mergeinit|
% Calculate the initial parameters for a merged mixture component to be
% substituted for two given components.
%
% *Syntax*
%
%   theta = D.mergeinit(idx1, idx2, theta)
%   theta = D.mergeinit(idx1, idx2, theta, options)
%   theta = D.mergeinit(idx1, idx2, theta, options, data)
%   [theta, store] = D.mergeinit(idx1, idx2, theta, options, data, store)
%
% *Description*
%
% |theta = D.mergeinit(idx1, idx2, theta)| returns the parameters for a new
% mixture resulted by merging the components indexed by |idx1| and |idx2|.
% The input |theta| is the parameters for the mixture before the merge.
% Initialization for the merged component is performed using the default
% methods for each parameter.
%
% |theta = D.mergeinit(idx1, idx2, theta, options)| uses
% |options.mergeInit| to initialize the parameters for the merged component
% as described <../estimation_options.html#10 here>.
%
% |theta = D.mergeinit(idx1, idx2, theta, options, data)| can be used if
% the initialization method requires estimation data.
%
% |[theta, store] = D.mergeinit(idx1, idx2, theta, options, data, store)|
% can be used for caching purposes as described <../caching.html here>.
%
% *Note:* This function is automatically called by the |merge| function and
% you don't usually need to call this function separately.
%

    D.mergeinit = @mergeinit;
    function [varargout] = mergeinit(varargin)
        [varargout{1:nargout}] = mixture_mergeinit(D, varargin{:}); 
    end

%% |merge|
% Merge two variable components into one
%
% *Syntax*
%
%   newD = D.merge(idx1, idx2)
%   [newD, newtheta] = D.merge(idx1, idx2, theta)
%   [newD, newtheta] = D.merge(idx1, idx2, theta, options)
%   [newD, newtheta] = D.merge(idx1, idx2, theta, options, data)
%   [newD, newtheta, idxMerged, idxMap] = D.merge(...)
%   [newD, newtheta, idxMerged, idxMap, store] = D.merge(idx1, idx2, theta, options, data, store)
%
% *Description*
%
% |newD = D.merge(idx1, idx2)| returns the new mixture, |newD|, resulted by
% merging the components indexed by |idx1| and |idx2| from |D|.
%
% |[newD, newtheta] = D.merge(idx1, idx2, theta)| also initializes the
% parameters for the merged components. The input |theta| is the parameter
% structure for the mixture |D| before the merge. The output |newtheta| is
% the parameter structure of the new mixture, |newD|, after the merge (with
% one less item than |theta|). The initialized parameter values for the
% merged component are stored at |newtheta.D{idx1}| and |theta.D{idx2}| is
% removed in |newtheta|. The component weights in |newtheta.p| are also
% updated accordingly, such that the merged component has the sum of the
% weights of the original components.
%
% |[newD, newtheta] = D.merge(idx1, idx2, theta, options)| uses
% |options.mergeInit| to initialize the parameters for the merged component
% as described <../estimation_options.html#10 here>.
%
% |[newD, newtheta] = D.merge(idx1, idx2, theta, options, data)| can be
% used if the initialization method requires estimation data.
%
% |[newD, newtheta, idxMerged, idxMap] = D.merge(...)| also returns the
% index of the merged component in |idxMerged|. The output |idxMap| is an
% index-map vector, with length |D.num()|, which maps the old indices of
% the components in |D| to their new indices in |newD|. You can use the map
% like this: |idx_new = idxMap(idx_old)|, where |idx_old| denotes an index
% for a component in |D| and |idx_new| denotes its new index in the
% components of |newD|.
%
% |[newD, newtheta, idxMerged, idxMap, store] = D.merge(idx1, idx2, theta,
% options, data, store)| can be used for caching purposes as described
% <../caching.html here>.
%
% *Example*
%
% The following code demonstrates a basic usage of merge-related functions:
%
%   % Construct a mixture distribution and generate random 
%   % parameters and data for demonstration
%   D = mixturefactory(mvnfactory(2), 3);
%   theta = D.randparam();
%   data = D.sample(theta, 1000);
%   % Find candidate components for merging
%   [idx_merge1, idx_merge2] = D.mergecandidates(theta, data);
%   % Merge the first candidate components
%   idx1 = idx_merge1(1);
%   idx2 = idx_merge2(1);
%   [newD, newtheta, idxMerged] = D.merge(idx1, idx2, theta);
%   % We may estimate the parameters merely for the merged component
%   newtheta = newD.estimatepartial(idxMerged, newtheta, data)
%

    D.merge = @merge;
    function [newD, newtheta, idxMerged, idxMap, store] = merge(idx1, idx2, theta, options, data, store)
        
        % sort given indices
        if idx1 > idx2
            temp = idx1;
            idx1 = idx2;
            idx2 = temp;
        end
        
        if nargout > 1
            if nargin < 6
                store = struct;
            end
            if nargin < 5
                data = [];
            end
            if nargin < 4
                options = mxe_options();
            else
                options = mxe_options(options);
            end
            [newtheta, store] = mergeinit(idx1, idx2, theta, options, data, store);
        end
        
        % calculate new indices
        if nargout > 2
            idxMerged = idx1;
        end
        
        % calculate the index map for original (variable) components.
        % can be used like:   idx_new = idxMap(idx_old)
        if nargout > 3
            idxMap = 1:num;
            idxMap(idx2) = idx1;
            idxMap(idxMap > idx2) = idxMap(idxMap > idx2) - 1;
        end
        
        % build new varD
        newVarD = varD;
        newVarD(idx2) = [];
        
        % build fixing struct
        Afixing = struct('fixing__', true);
        Afixing.varD = newVarD;
        Afixing.fixedD = fixedD;
        Afixing.fixedtheta = fixedtheta;
        if cacheValid
            Afixing.hXcache = hXcache; % cache is still valid
        end
        
        % construct the new distribution
        newD = mixturefactory(Afixing);
    end

%% |splitcandidates|
% Find split candidates
%
% *Syntax*
%
%   idx = D.splitcandidates(theta, data)
%   idx = D.splitcandidates(theta, data, options)
%   idx = D.splitcandidates(theta, data, options, n)
%
% *Description*
%
% |idx = D.splitcandidates(theta, data)| returns the vector |idx|
% containing indices of the components of |D|, sorted based on the default
% split criterion (see <../estimation_options.html#10 Split-and-merge
% options>). |idx| has a length equal to the number of components
% (|D.num()|). The first element in |idx| refers to the best candidate
% component for splitting, and the last element refers to the one having
% the worse value of the split criterion.
%
% |idx = D.splitcandidates(theta, data, options)| uses
% |options.sm.splitCriterion| as the split criterion for sorting the
% components. see <../estimation_options.html#10 Split-and-merge options>
% for more details.
%
% |idx = D.splitcandidates(theta, data, options, n)| where |n| is a
% positive integer less than or equal to |D.num()|, returns only the |n|
% best split candidates (|idx| will have a length of |n|).
%
% *Example*
%
% The following code demonstrates a basic usage of split-related functions:
%
%   % Construct a mixture distribution and generate random 
%   % parameters and data for demonstration
%   D = mixturefactory(mvnfactory(2), 3);
%   theta = D.randparam();
%   data = D.sample(theta, 1000);
%   % Find candidate components for splitting
%   idx_split = D.splitcandidates(theta, data);
%   % Split the first candidate component
%   idx = idx_split(1);
%   [newD, newtheta, idxSplitted] = D.split(idx, theta);
%   % We may estimate the parameters merely for the splitted components
%   newtheta = newD.estimatepartial(idxSplitted, newtheta, data)
%

    D.splitcandidates = @splitcandidates;
    function [varargout] = splitcandidates(varargin)
        [varargout{1:nargout}] = mixture_splitcandidates(D, varargin{:}); 
    end

%% |mergecandidates|
% Find merge candidates
%
% *Syntax*
%
%   [idx1, idx2] = D.mergecandidates(theta, data)
%   [idx1, idx2] = D.mergecandidates(theta, data, options)
%   [idx1, idx2] = D.mergecandidates(theta, data, options, n)
%
% *Description*
%
% |[idx1, idx2] = D.mergecandidates(theta, data)| returns the vectors
% |idx1| and |idx2| containing index pairs of the components of |D|, sorted
% based on the default merge criterion (see <../estimation_options.html#10
% Split-and-merge options>). |idx1| and |idx2| both have a length equal to
% the number of components (|D.num()|). The first elements in |idx1| and
% |idx2| refer to the best pair of candidate components for merging, and
% the last pair of elements refer to those having the worse value of the
% merge criterion.
%
% |[idx1, idx2] = D.mergecandidates(theta, data, options)| uses
% |options.sm.mergeCriterion| as the merge criterion for sorting the
% components. see <../estimation_options.html#10 Split-and-merge options>
% for more details.
%
% |[idx1, idx2] = D.mergecandidates(theta, data, options, n)| where |n| is
% a positive integer less than or equal to |D.num()|, returns only the |n|
% best merge candidates (|idx1| and |idx2| will each have a length of |n|).
%
% *Example*
%
% The following code demonstrates a basic usage of merge-related functions:
%
%   % Construct a mixture distribution and generate random 
%   % parameters and data for demonstration
%   D = mixturefactory(mvnfactory(2), 3);
%   theta = D.randparam();
%   data = D.sample(theta, 1000);
%   % Find candidate components for merging
%   [idx_merge1, idx_merge2] = D.mergecandidates(theta, data);
%   % Merge the first candidate components
%   idx1 = idx_merge1(1);
%   idx2 = idx_merge2(1);
%   [newD, newtheta, idxMerged] = D.merge(idx1, idx2, theta);
%   % We may estimate the parameters merely for the merged component
%   newtheta = newD.estimatepartial(idxMerged, newtheta, data)
%

    D.mergecandidates = @mergecandidates;
    function [varargout] = mergecandidates(varargin)
        [varargout{1:nargout}] = mixture_mergecandidates(D, varargin{:}); 
    end

%% |dim|
% See <doc_distribution_common.html#3 distribution structure common members>.

    D.dim = @dim; % parameter space dimensions
    function dim = dim()
        dim = nump - 1; % component weights
        for k = 1:num
            dim = dim + Components{k}.dim(); % component parameters
        end
    end

%% |datadim|
% See <doc_distribution_common.html#4 distribution structure common members>.

    D.datadim = @() Components{1}.datadim(); % data space dimensions

%%

    function store = weighting_intermediate_params(theta, data, store, datapatchsize)
    % calculate intermediate parameters for weighting
        
        if ~isfield(store, 'componentStores')
            store.componentStores = cell(num,1);
        end
        
        if nargin < 4
            datapatchsize = inf;
        end
        
        if ~isfield(store, 'hX') || ~isfield(store, 'llik')
            data = mxe_readdata(data, false);
            index = data.index;
            
            % add fixed lls to hX before calculating total ll
            if numfixed > 0
                if cacheValid
                    hX_fixed = hXcache(:,index);
                else
                    hX_fixed = calc_cache(data, datapatchsize);
                end
            end
            
            n = data.size;
            data.weight = [];
            
            % Initialize different variables
            hX = zeros(num, n);
            
            if isinf(datapatchsize)
                datapatchsize = n;
            end
            
            numk2 = ceil(n / datapatchsize);
            
           % Calculate the log-likelihood of data goes to different clusters
            for k2 = 1:numk2
                
                bind = 1 + (k2-1)*datapatchsize;
                eind = min(k2*datapatchsize, data.size);
                data.index = index(bind:eind); % done with patching
                for k = 1:num
                    if ~isstruct(store.componentStores{k}) || numk2 > 1
                        store.componentStores{k} = struct;
                    end
                    [llvec, store.componentStores{k}] = ...
                        Components{k}.llvec(...
                        theta.D{k}, data, store.componentStores{k});
                    hX(k,bind:eind) = log(theta.p(k)) + llvec;
                end
                
            end
            
           
            
            % add fixed lls to hX before calculating total ll
            if numfixed > 0
                hX = vertcat(hX, log(theta.p(num+1)) + hX_fixed);
            end
            
            % hX is a (num-by-n) matrix where each column of hX contains
            % log(p_k*D_k) for a data point for every component (1<k<num).
            store.hX = hX;
            % Based on log-likelihoods calculate the total log-likelihood and
            % weights
            store.llik = logsumexp(hX, 1);
        end
    end

%% |weighting|
% Calculate the weighting (posterior probability) of each mixture component
% for each data point
%
% *Syntax*
%
%   component_weights = D.weighting(theta, data)
%   [component_weights, store] = D.weighting(theta, data, store)
%
% *Description*
%
% |component_weights = D.weighting(theta, data)| returns a |K-by-N| matrix
% where |K| is the number of mixture components (|D.num()|) and |N| is the
% number of data points. Each column of the matrix contains the posterior
% probabilities of each mixture component for the corresponding data point.
% The probabilities are normalized such that each column sums up to unity.
% |theta| is the parameter structure for the mixture, and |data| is the
% data input.
%
% |[component_weights, store] = D.weighting(theta, data, store)| can be
% used for caching purposes as described <../caching.html here>.
%
% For information about the input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *Example*
%
%   D = mixturefactory(mvnfactory(1), 2);
%   theta.D{1} = struct('mu',0, 'sigma',1);
%   theta.D{2} = struct('mu',2, 'sigma',1);
%   theta.p = [0.5 0.5];
%   data = [0 1 2];
%   component_weights = D.weighting(theta, data)
%
%   component_weights =
%     0.8808    0.5000    0.1192
%     0.1192    0.5000    0.8808
%

    D.weighting = @weighting;
    function [component_weights, store] = weighting(theta, data, store, datapatchsize)
        
        if nargin < 3
            store = struct;
        end
        
        if nargin < 4
            datapatchsize = inf;
        end
        
        store = weighting_intermediate_params(theta, data, store, datapatchsize);
        hX = store.hX;
        llik = store.llik;
        
        weight = mxe_readweight(data);
            
        % undo logarithm and normalize hX such that each column sums up to 1.
        component_weights = exp( bsxfun(@minus, hX, llik) );
        
        if ~isempty(weight)
            component_weights = bsxfun(@times, component_weights, weight);
        end
        low_memory = true;
        if low_memory
            store = rmfield(store,'hX');
            %store = rmfield(store,'llik');
        end
    end

%% |ll|
% See <doc_distribution_common.html#5 distribution structure common members>.

    D.ll = @ll;
    function [ll, store] = ll(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        store = weighting_intermediate_params(theta, data, store);
        llik = store.llik;
        
        weight = mxe_readweight(data);
      
        % sum the log-likelihoods
        if ~isempty(weight)
            ll = sum(weight .* llik);
        else
            ll = sum(llik);
        end
        
        if low_memory
            store = rmfield(store,'hX');
            store = rmfield(store,'llik');
        end
        
    end

%% |llvec|
% See <doc_distribution_common.html#6 distribution structure common members>.

    D.llvec = @llvec;
    function [llvec, store] = llvec(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        store = weighting_intermediate_params(theta, data, store);
        llik = store.llik;
        
        weight = mxe_readweight(data);
      
        % sum the log-likelihoods
        if ~isempty(weight)
            llvec = weight .* llik;
        else
            llvec = llik;
        end
        
        if low_memory
            store = rmfield(store,'hX');
            store = rmfield(store,'llik');
        end
        
    end

%% |llgrad|
% See <doc_distribution_common.html#7 distribution structure common members>.

    D.llgrad = @llgrad;
    function [dll, store] = llgrad(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        [component_weights, store] = weighting(theta, data, store);

        data = mxe_readdata(data, false);
        index = data.index;
        datamat = data.data;

        % Given the weighting calculate the gradient
        dll.D = cell(num, 1);
        for k = 1:num
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            [dll.D{k}, store.componentStores{k}] = ...
                Components{k}.llgrad(theta.D{k}, ...
                struct('data',datamat, 'weight',component_weights(k,:), 'index',index), ...
                store.componentStores{k});
        end
        dll.p = squeeze(sum(component_weights, 2)./theta.p); %TODO: this uses the fact that component_weights and theta.p have the same number of rows, accidentally

    end

%% |llgraddata|
% See <doc_distribution_common.html#8 distribution structure common members>.

    D.llgraddata = @llgraddata;
    function [dld, store] = llgraddata(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        [component_weights, store] = weighting(theta, data, store);
        data = mxe_readdata(data, false);
        index = data.index;
        datamat = data.data;

        % Given the weighting calculate the gradient
        dld = 0;
        for k = 1:num
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            [dltemp, store.componentStores{k}] = ...
                Components{k}.llgraddata(theta.D{k}, ...
                struct('data',datamat, 'weight',component_weights(k,:), 'index',index), ...
                store.componentStores{k});
            dld = dld + dltemp;
        end
        
    end

%% |cdf|
% See <doc_distribution_common.html#9 distribution structure common members>.

    D.cdf = @cdf;
    function y = cdf(theta, data)
        data = mxe_readdata(data);
        N = data.size;

        y = zeros(1, N);
        fulltheta = fullparam(theta);
        for k = 1:numtotal
            y = y + fulltheta.p(k) .* Components{k}.cdf(fulltheta.D{k}, data);
        end
    end

%% |pdf|
% See <doc_distribution_common.html#10 distribution structure common members>.

%% |sample|
% See <doc_distribution_common.html#11 distribution structure common members>.
%
% You can also get the component labels for each data point using the
% following syntax:
%
%   [data, label] = D.sample(...)
%

    D.sample = @sample;
    function [data, label] = sample(theta, n)
        
        if nargin < 2, n = 1; end
        
        q = D.datadim();
        data = zeros(q, n);
        label = zeros(1, n);
        ind = randsample(1:numtotal,n,true,theta.p);
        fulltheta = fullparam(theta);
        for k = 1:numtotal
            data(:, ind == k) = Components{k}.sample(fulltheta.D{k}, sum(ind == k));
            label(1, ind == k) = k;
        end
    end

%% |randparam|
% See <doc_distribution_common.html#12 distribution structure common members>.

%% |init|
% See <doc_distribution_common.html#13 distribution structure common members>.

    D.init = @init;
    function theta = init(data, varargin)

        data = mxe_readdata(data);
        label = ceil(rand(1, data.size) * num);
        data = data.data;
        theta.D = cell(num,1);
        for k = 1:num
            theta.D{k} = Components{k}.init(data(:,label == k), varargin{:});
        end
        theta.p = ones(nump, 1) ./ nump;
    end

%% |estimatedefault|
% Default estimation function for mixture distribution. This function
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
% *Example*
%
%   % create a mixture of two Gaussian distributions
%   D = mixturefactory(mvnfactory(1), 2);
%   % generate 1000 random data points
%   data = [randn(1,300).*2-5, randn(1,700).*3+10];
%   % use an appropriate theta0
%   options.theta0 = D.init(data);
%   % estimate mixture parameters using the EM algorithm
%   theta = D.estimatedefault(data, options)
%

    D.estimatedefault = @estimatedefault;
    function [varargout] = estimatedefault(varargin)
        [varargout{1:nargout}] = mixture_estimatedefault(D, varargin{:}); 
    end

%% |estimatepartial|
% Estimate parameters for a subset of components, while fixing the others
%
% *Syntax*
%
%   newtheta = D.estimatepartial(idx, theta, data)
%   newtheta = D.estimatepartial(idx, theta, data, options)
%   [newtheta, D] = D.estimatepartial(...)
%   [newtheta, D, info] = D.estimatepartial(...)
%   [newtheta, D, info, options] = D.estimatepartial(...)
%
% *Description*
%
% |newtheta = D.estimatepartial(idx, theta, data)| estimates the parameters
% for the components indexed by the vector |idx| while preserving the
% parameter values for the other mixture components. |theta| contains the
% current values of parameters for all mixture components. The parameter
% values for those components that are indexed by |idx| are used as the
% initial point for their estimation, and the other parameters remain
% unchanged in the returned parameter structure.
%
% |newtheta = D.estimatepartial(idx, theta, data, options)| utilizes
% applicable options from the |options| structure in the estimation
% procedure.
%
% |[newtheta, D] = D.estimatepartial(...)| also returns |D|, the
% distribution structure for which |theta| is applicable. (This is the same
% as the distribution structure |D| from which you called |estimate|, and
% so it should not normally be used. The purpose of including it in the
% output is to maintain compatibility with other estimation functions).
%
% |[newtheta, D, info] = D.estimatepartial(...)| also returns |info|, a
% structure array containing information about successive iterations
% performed by iterative estimation functions.
%
% |[newtheta, D, info, options] = D.estimatepartial(...)| also returns the
% effective |options| used, so you can see what default values the function
% used on top of the options you possibly specified.
%
% For information about the input |theta| and output |newtheta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>. You may also want to read about
% <../estimation_options.html |options|> or
% <../estimation_statistics_structure.html |info|> arguments.
%
% *Available Options*
%
% This function supports all the options described in
% <../estimation_options.html estimation options>.
%
% *Note:* |options.theta0| will be ignored since the input |theta| contains
% the initial point for estimation.
%
% *Returned |info| fields*
%
% The fields present in the returned |info| structure array, depend on the
% solver used (|options.solver|). When a Manopt solver is specified, the
% |info| returned by the Manopt solver is returned directly. For the 'default'
% solver see the documentation of the 'estimatedefault' function for the
% specific distribution. You can read more at our documentation on
% <../estimation_statistics_structure.html estimation statistics
% structure>.
%
% *Example*
%
% The following code uses |estimatepartial| to estimate the parameters of
% splitted components after a split:
%
%   % Construct a mixture distribution and generate random 
%   % parameters and data for demonstration
%   D = mixturefactory(mvnfactory(2), 3);
%   theta = D.randparam();
%   data = D.sample(theta, 1000);
%   % Find candidate components for splitting
%   idx_split = D.splitcandidates(theta, data);
%   % Split the first candidate component
%   idx = idx_split(1);
%   [newD, newtheta, idxSplitted] = D.split(idx, theta);
%   % Estimate the parameters merely for the splitted components
%   newtheta = newD.estimatepartial(idxSplitted, newtheta, data)
%

    D.estimatepartial = @estimatepartial;
    function [newtheta, newD, info, options] = estimatepartial(idx, theta, data, options)

        if nargin < 4
            options = mxe_options();
        else
            options = mxe_options(options);
        end

        % make the other components fixed
        invidx = invertindex(idx);
        [newD, theta0, idxfixed] = fixate(invidx, theta, data, options.datapatchsize);
        
        % use the given theta as the initial point for estimation
        options.theta0 = theta0;
        
        % estimate the resulting partial mixture parameters
        [newtheta, newD, info, options] = newD.estimate(data, options);
        
        % bring the fixed components back
        [newD, newtheta] = newD.unfix(idxfixed, invidx, newtheta);
    end

%% |penalizerparam|
% See <doc_distribution_common.html#15 distribution structure common members>.
%
% *Penalizer Info*
%
% The default penalizer for the mixture distribution is the sum of the
% default penalizers of its components.
%

    D.penalizerparam = @penalizerparam;
    function penalizer_theta = penalizerparam(data)
        penalizer_theta.phi = ones(nump,1);
        for k = 1:num
            penalizer_theta.D{k} = Components{k}.penalizerparam(data);
        end
    end

%% |penalizercost|
% See <doc_distribution_common.html#16 distribution structure common members>.

    D.penalizercost = @penalizercost;
    function [costP, store] = penalizercost(theta, penalizer_theta, store)
        
        if nargin < 3
            store = struct;
        end
        
        if ~isfield(store, 'componentStores')
            store.componentStores = cell(num,1);
        end
        
        costP = sum((penalizer_theta.phi-1).*log(theta.p));
        for k = 1:num
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            [cost, store.componentStores{k}] = ...
                Components{k}.penalizercost(theta.D{k}, penalizer_theta.D{k}, ...
                store.componentStores{k});
            costP = costP + cost;
        end
    end

%% |penalizergrad|
% See <doc_distribution_common.html#17 distribution structure common members>.

    D.penalizergrad = @penalizergrad;
    function [gradP, store] = penalizergrad(theta, penalizer_theta, store)
        
        if nargin < 3
            store = struct;
        end
        
        if ~isfield(store, 'componentStores')
            store.componentStores = cell(num,1);
        end
        
        gradP.p = (penalizer_theta.phi-1)./theta.p;
        for k = 1:num
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            [gradP.D{k}, store.componentStores{k}] = ...
                Components{k}.penalizergrad(theta.D{k}, penalizer_theta.D{k}, ...
                store.componentStores{k});
        end
    end

%% |regcost|
% Regularizer cost
%
% *Syntax*
%
%   reg = D.regcost(theta, data)
%   [reg, store] = D.regcost(theta, data, store)
%
% *Description*
%
% |reg = D.regcost(theta, data)| returns the cost penalty for the
% model-free regularization method (Hosseini, 2012). This regularization
% penalizes components that take few data points.
%
% |[reg, store] = D.regcost(theta, data, store)| can be used for caching
% purposes as described <../caching.html here>.
%
% The output of this function is used as a regularizer for the cost
% function during parameter estimation when |options.regularize| is turned
% on.
%
% For information about the input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *References*
%
% # R. Hosseini, “Natural Image Modelling using Mixture Models with
% compression as an application,” Berlin, Technische Universtit?t Berlin,
% Diss., 2012, 2012.
%

    D.regcost = @regcost;
    function [reg, store] = regcost(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        [component_weights, store] = weighting(theta, data, store);

        % First compute the number of data that go to different clusters
        if ~isfield(store, 'hXsum')
            store.hXsum = sum(component_weights, 2);
        end
        
        % Compute the cost with respect to component weights
        regparam = dim();
        regparam = max(regparam , sum(store.hXsum)/(10*num));
        reg = BayesFactor(store.hXsum, regparam);
        reg2 = HAICmc(store.hXsum, regparam);
        reg = reg - reg2;
    end

%% |reggrad|
% Regularizer gradient with respect to parameters
%
% *Syntax*
%
%   dll = D.reggrad(theta, data)
%   [dll, store] = D.reggrad(theta, data, store)
%
% *Description*
%
% |dll = D.reggrad(theta, data)| returns the (Euclidean) gradient penalty
% for the model-free regularization method (Hosseini, 2012). This
% regularization penalizes components that take few data points.
%
% |[dll, store] = D.reggrad(theta, data, store)| can be used for caching
% purposes as described <../caching.html here>.
%
% The output of this function is used as a regularizer for the
% cost-gradient function during parameter estimation when
% |options.regularize| is turned on.
%
% For information about the input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *References*
%
% # R. Hosseini, “Natural Image Modelling using Mixture Models with
% compression as an application,” Berlin, Technische Universtit?t Berlin,
% Diss., 2012, 2012.
%

    D.reggrad = @reggrad;
    function [dll, store] = reggrad(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        [component_weights, store] = weighting(theta, data, store);

        % First compute the number of data that go to different clusters
        if ~isfield(store, 'hXsum')
            store.hXsum = sum(component_weights, 2);
        end
        
        % Compute the gradient with respect to component weights
        regparam = dim();
        regparam = max(regparam , sum(store.hXsum)/(10*num));
        [unused, dreg] = BayesFactor(store.hXsum, regparam); %#ok<ASGLU>
        [unused, dreg2] = HAICmc(store.hXsum, regparam); %#ok<ASGLU>
        dreg = dreg - dreg2;
        % based on regularization gradient update the component weighting
        weight_mod = bsxfun(@times, dreg, component_weights);
        common_part = col_sum(weight_mod);
        
        component_weights = component_weights+ weight_mod - ...
            bsxfun(@times, component_weights, common_part);
        
        data = mxe_readdata(data, false);
        index = data.index;
        datamat = data.data;

        % Given the weighting calculate the gradient
        dll.D = cell(num, 1);
        for k = 1:num
            if ~isstruct(store.componentStores{k})
                store.componentStores{k} = struct;
            end
            [dll.D{k}, store.componentStores{k}] = ...
                Components{k}.llgrad(theta.D{k}, ...
                struct('data',datamat, 'weight',component_weights(k,:), 'index',index), ...
                store.componentStores{k});
        end
        dll.p = squeeze(sum(component_weights, 2)./theta.p);
    end

%% |sumparam|
% See <doc_distribution_common.html#18 distribution structure common members>.

    D.sumparam = @sumparam;
    function theta = sumparam(theta1, theta2)
        theta.D = cell(num, 1);
        for k = 1:num
            theta.D{k} = Components{k}.sumparam(theta1.D{k}, theta2.D{k});
        end
        theta.p = theta1.p + theta2.p; %because it is used in gradient
        %theta.p = theta.p ./ sum(theta.p); %incorrect
    end

%% |scaleparam|
% See <doc_distribution_common.html#19 distribution structure common members>.

    D.scaleparam = @scaleparam;
    function theta = scaleparam(scalar, theta)
        for k = 1:num
            theta.D{k} = Components{k}.scaleparam(scalar, theta.D{k});
        end
        theta.p = scalar * theta.p; %because it is used in gradient
%         theta.p = theta.p; incorrect
    end

%% |sumgrad|
% See <doc_distribution_common.html#20 distribution structure common members>.

%% |scalegrad|
% See <doc_distribution_common.html#21 distribution structure common members>.

%% |entropy|
% See <doc_distribution_common.html#22 distribution structure common members>.

    D.entropy = @entropy;
    function h = entropy(theta)
        h = zeros(1,numtotal);
        fulltheta = fullparam(theta);
        for k = 1:numtotal
            h(k) = Components{k}.entropy(fulltheta.D{k});
        end
    end

%% |kl|
% Calculate Kullback–Leibler divergence between mixture components
%
% *Syntax*
%
%   kl = D.kl(theta)
%
% *Description*
%
% |kl = D.kl(theta)| returns a symmetric |K-by-K| matrix where |K| is the
% number of mixture components (|D.num()|). The element at (i,j) in this
% matrix contains the KL divergence between the i'th and j'th components.
% |theta| is the parameter structure for the mixture.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>.
%

    D.kl = @kl;
    function kl = kl(theta)
        kl = zeros(numtotal);
        fulltheta = fullparam(theta);
        for k1 = 1:numtotal
            for k2 = 1:k1-1
                kl(k1,k2) = Components{1}.kl(fulltheta.D{k1}, fulltheta.D{k2});
                kl(k2,k1) = kl(k1,k2);
            end
        end        
    end

%% |AICc|
% See <doc_distribution_common.html#24 distribution structure common members>.

%% |BIC|
% See <doc_distribution_common.html#25 distribution structure common members>.

%%

    function [hYsum, store] = sumcluster(theta, data, store)
    % The data that go to a specific cluster
    
        if nargin < 3
            store = struct;
        end
        
        [hX, store] = weighting(theta, data, store);
        hYsum = sum(hX, 2);
    end

%% |MML|
% Calculate minimum message length information criterion
%
% *Syntax*
%
%   mml = D.MML(theta, data)
%   [mml, gradMML] = D.MML(theta, data)
%
% *Description*
%
% |mml = D.MML(theta, data)| calculates the minimum message length (MML)
% information criterion of the distribution |D| with parameters |theta| for
% the given |data|.
%
% |[mml, gradMML] = D.MML(theta, data)| also returns |gradMML|, the
% gradient of MML with respect to component posterior probabilities.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%

    D.MML = @MML;
    function [mml, gradMML] = MML(theta, data)
    
        regparam = dim();
        hYsum = sumcluster(theta, data);
        if nargout == 1
            mml = BayesFactor(hYsum, regparam);
        else
            [mml, gradMML] = BayesFactor(hYsum, regparam);
        end
        n = sum(hYsum);
        mml = mml + regparam/2 * log(n) + regparam/2*(1-log(12)); 
    end
    function [cost, gradweight] = BayesFactor(hYn, regparam)
    %
    %  [cost,gradweight] = BICm(hYn,regparam)
    %
    % It calculates the bayes factor for the mixture model
    %
    %  Inputs:
    %    hYn: The number of data goes to different clusters
    %    regparam  : total number of parameters in the mixture model
    %
    %  Outputs:
    %    cost: cost function
    %    gradweight: gradient of the cost function with respect to hYn
    %

        num_com = length(hYn);
        regparam = regparam /  num_com;
        cost = regparam/2 * sum(sum(log(hYn))) - regparam/2 * log(2*pi);
        if nargout > 1
            gradweight = regparam/2 ./ hYn;
        end
    end

%% |BICm|
% Calculate modified Bayesian information criterion
%
% *Syntax*
%
%   bic = D.BICm(theta, data)
%
% *Description*
%
% |bic = D.BICm(theta, data)| calculates the modified Bayesian information
% criterion (Hosseini, 2012) of the distribution |D| with parameters
% |theta| for the given |data|.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *References*
%
% # R. Hosseini, “Natural Image Modelling using Mixture Models with
% compression as an application,” Berlin, Technische Universtit?t Berlin,
% Diss., 2012, 2012.
%

    D.BICm = @BICm;
    function bic = BICm(theta, data)
    %
    %  [bic, gradBic] = BICm(theta,dat)
    %
    % The function calculated  BIC information criteria for the model 
    %
    %  Inputs:
    %    theta  : theta is mcgsm distribution class
    %    Data: A NISDET class variable for Data
    %
    %  Outputs:
    %    aic: the total criteria value
    %    gradBic: A vector that contains the gradient with respect to data
    %    that go to different clusters
    %

        regparam = dim();
        hYsum = sumcluster(theta, data);
        bic = BayesFactor(hYsum, regparam);
        n = sum(hYsum);
        bic = bic + regparam/2 * log(n) - regparam*log(2*pi); 
    end

%% |AICmc|
% Calculate modified corrected Akaike information criterion
%
% *Syntax*
%
%   aicc = D.AICmc(theta, data)
%   [aicc, gradAicc] = D.AICmc(theta, data)
%
% *Description*
%
% |aicc = D.AICmc(theta, data)| calculates the modified corrected Akaike
% information criterion (Hosseini, 2012) of the distribution |D| with
% parameters |theta| for the given |data|.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *References*
%
% # R. Hosseini, “Natural Image Modelling using Mixture Models with
% compression as an application,” Berlin, Technische Universtit?t Berlin,
% Diss., 2012, 2012.
%

    D.AICmc = @AICmc;
    function [aicc, gradAicc] = AICmc(theta, data)
    %
    %  [aicc,gradAicc]=icAICc(theta,dat)
    %
    % The function calculates AICc information criteria for the model
    %
    %  Inputs:
    %    theta  : theta is mcgsm distribution class
    %    Data: A NISDET class variable for Data
    %
    %  Outputs:
    %    aic: the total criteria value
    %    gradAicc: A vector that contains the gradient with respect to data
    %    that go to different clusters
    %

        regparam = dim();
        hYsum = sumcluster(theta, data);
        if nargout == 1
            aicc = HAICmc(hYsum, regparam);
        else
            [aicc, gradAicc] = HAICmc(hYsum, regparam);
        end
    end
    function [cost, gradweight] = HAICmc(hYn, regparam)
    %
    %  [cost,gradweight] = AICmc(hYn,regparam,predparam)
    %
    % It calculates the Huvrich et al. corrected AIC for mixture model
    %
    %  Inputs:
    %    hYn: The number of data goes to different clusters
    %    regparam  : total number of parameters in the mixture model
    %
    %  Outputs:
    %    cost: cost function
    %    gradweight: gradient of the cost function with respect to hYn
    %

        num_com = length(hYn);

        param = regparam / num_com;
        if nargout == 1
            y = func(hYn - param - 1);
            cost = param * sum(sum(hYn ./ y));
            if any(hYn < param + 1)
                cost = 1e15;
            end
                
        else
            [y, dy] = func(hYn - param - 1);
            cost = param * sum(sum(hYn ./ y));
            % the gradient
            gradweight = param * (y - dy.*hYn) ./ y.^2;
            if any(hYn < param + 1)
                cost = 1e15;
                gradweight = ones(size(gradweight)) * 1e15;
            end
        end
    end
    function [y, dy] = func(x)
        a = 5;
        y = 1/a * log(1 + exp(a * x));
        y(x > 100) = x(x > 100);
        y(x < -2) = 1/a * exp(a * x(x < -2));
        if nargout > 1
            dy = 1 ./ (1 + exp(-a * x));
            dy(x > 100) = 1;
        end
        dy(x < -2) = exp(a * x(x < -2));
    end

%% |display|
% See <doc_distribution_common.html#26 distribution structure common members>.

    D.display = @display;
    function str = display(theta)
        str = '';
        for k = 1:num
            str = [str, sprintf('\nvarD(%d): %s / weight: %g\n', k, varD{k}.name(), theta.p(k))]; %#ok<AGROW>
            str = [str, varD{k}.display(theta.D{k})]; %#ok<AGROW>
        end
        if numfixed > 0
            w = theta.p(num+1);
            for k = 1:numfixed
                str = [str, sprintf('\nfixedD(%d): %s / weight: %g\n', k, fixedD{k}.name(), w.*fixedtheta.p(k))]; %#ok<AGROW>
                str = [str, fixedD{k}.display(fixedtheta.D{k})]; %#ok<AGROW>
            end
        end
        
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
        [varargout{1:nargout}] = mixture_visualize(D, varargin{:});
    end

%% |gaussianize|
%
% *Syntax*
%
%   y = gaussianize(theta, data)
%

    D.gaussianize = @gaussianize;
    function y = gaussianize(theta, data)
        y = cdf(theta, data);
        for k = 1: size(y,1)
            y(k,:) = -sqrt(2).*erfcinv(2*y(k,:));
        end
    end


%%

    D = mxe_addsharedfields(D);
end