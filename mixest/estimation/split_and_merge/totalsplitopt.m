%% |totalsplitopt|
% Total-split algorithm for estimation of mixture model structure and
% parameters, with generalized internal estimations.
%
% *Syntax*
%
%   [theta, D] = totalsplitopt(data, target_num)
%   [theta, D] = totalsplitopt(data, target_num, options)
%   [theta, D, info] = totalsplitopt(...)
%   [theta, D, info, options] = totalsplitopt(...)
%
% *Description*
%
% This function implements the total-split (split-then-merge) algorithm
% proposed by Sra et al. (2015), though here we use generalized
% optimization methods.
%
% |[theta, D] = totalsplit(data, target_num)| returns the estimated
% parameters and the mixture distribution |D| with |target_num| components,
% fitted to |data|.
%
% |[theta, D] = totalsplit(data, target_num, options)| utilizes applicable
% options from the |options| structure in the estimation procedure.
%
% |[theta, D, info] = totalsplit(...)| also returns |info|, a structure
% array containing information about successive iterations performed by
% iterative estimation functions.
%
% |[theta, D, info, options] = totalsplit(...)| also returns the effective
% |options| used, so you can see what default values the function used on
% top of the options you possibly specified.
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
% This function supports all the options described in
% <../estimation_options.html estimation options>. This function accepts
% the following additional fields in |options.sm|:
%
% * *|tolCostDiff|* (default |0|) : Minimum decrease in cost to accept a
% split.
% * *|numMax|* (default |2*target_num|) : Maximum number of mixture components.
% * *|componentD|* (default MVN distribution) : distribution structure
% defining the mixture component distribution type.
%
% |options.inner| can be used to set specific options for the inner
% estimations. To set options only for the partial or full inner
% estimations use |options.inner.partial| or |options.inner.full|
% respectively.
%
% *Returned |info| fields*
%
% The fields present in the returned |info| structure array, depend on the
% solver used (|options.solver|). When a Manopt solver is specified, the
% |info| returned by the Manopt solver is returned directly. For the 'default'
% solver see the documentation of the 'estimatedefault' function for the
% mixture distribution. You can read more at our documentation on
% <../estimation_statistics_structure.html estimation statistics
% structure>.
%
% *Example*
%
%   % generate 1000 random data points
%   data = randn(1,1000) .* 2 + 1;
%   % set some options
%   options.solver = 'cg';
%   options.verbosity = 2;
%   options.sm.numMax = 7;
%   options.inner.partial.maxiter = 10;
%   % fit mixture model to data
%   [theta, D] = totalsplit(data, 5, options)
%
% *References*
%
% # S. Sra, R. Hosseini, L. Theis, and M. Bethge, “Data modeling with the
% elliptical gamma distribution,” in Proceedings of the Eighteenth
% International Conference on Artificial Intelligence and Statistics, 2015,
% pp. 903–911.
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

function [theta, D, info, options] = totalsplit(data, target_num, options)
% This function implements the Total-Split method 
%
% This procedure applies splitting approach (partial splitting to increase
% the speed). Then it applies partial merging and once it reaches the
% suitable number of mixtures it stops. (If model selection is performed,
% then after each merge and optimization is applied) -> Since finding a
% mixture to split is hard (this is a good and fast strategy).
%
% During splitting, it can detect where singularity happens (using AICc)
% and avoids splitting that mixture further.
%

    extra_options = struct(...
        'sm', struct(...
            'costtype', 1, ...
            'tolcostdiff', 0, ...
            'nummax', ceil(1.6*target_num), ...
            'componentd', [], ... default is mvnfactory(datadim)
            'mixture', [], ...
            'savefolder', 'result'...
        ));

    if nargin < 3
        options = mxe_options([], extra_options);
    else
        options = mxe_options(options, extra_options);
    end
    
    if ~isempty(options.sm.mixture)
        D = options.sm.mixture; 
    elseif ~isempty(options.sm.componentd)
        ComponentD = options.sm.componentd; 
        D = mixturefactory(ComponentD, 1);
    else
        ds = mxe_readdata(data, false);
        ComponentD = mvnfactory(ds.dim);
        D = mixturefactory(ComponentD, 1);
        clear ds
    end
    
    
    % visualization init
    if options.visualization.enabled
        if ~isempty(options.visualization.visobject)
            vis = options.visualization.visobject;
        else
            vis = mxe_visualization(data, options.visualization);
            options.visualization.visobject = vis;
        end
    end

    % get the inner options for full and partial estimations
    opt_full = mxe_inneroptions(options, [], 'full');
    opt_partial = mxe_inneroptions(options, [], 'partial');
    
    % 1. start with a mixture of one distribution
    if options.verbosity >= 1
        fprintf('\nStarting with one component...\n')
    end
    
    [theta, D, info] = D.estimate(data, opt_full);
    
    % initializations
    cost = calc_cost(D, theta, data, options.sm.costtype, options);
    impvec = Inf; % vector of improvements achieved by the split that has produced each mixture component
    
    % First Loop : Splitting part
    while true
        
        % check if visualization figure is closed
        if options.visualization.enabled && options.visualization.stoponclose
            if vis.closed()
                return
            end
        end
        
        % Determine which mixture components should be splitted
        [unused, idx] = max(impvec); %#ok<ASGLU>
        
        % based on the index split the component
        if options.verbosity >= 1
            fprintf('\nSplitting component %d from %d mixture components...\n', idx, D.num())
            if options.debug
                fprintf('Improvement vector: %s\n', mat2str(impvec, 4))
%                 pause
            end
        end
        [newD, newtheta, idxSplitted, idxMap] = D.split(idx, theta, options, data);
        
        % Estimating the parameters partially for the newly splitted components
        opt = mxe_adjustoptions(opt_partial, info(end));
        opt.previnfo = info;
        [newtheta, newD, info] = newD.estimatepartial(idxSplitted, newtheta, data, opt);
        
        if newD.num() > 2
            % After partial estimation do also one round of full estimation
            opt_full.theta0 = newtheta;
            [newtheta, newD] = newD.estimate(data, opt_full);
        end
        
        namesave = [options.sm.savefolder '/tsplit' num2str(newD.num())];
        save(namesave, 'newtheta', 'newD');
        
        % calculate the cost improvement
        newcost = calc_cost(newD, newtheta, data, options.sm.costtype, options);
        costdiff = cost - newcost;
        
        if costdiff < options.sm.tolcostdiff
            % singularity detected: Reject the splitting and put end-node -Inf
            if options.verbosity >= 1
                fprintf('\nRejected split of component %d with cost improvement %g\n', idx, costdiff)
            end
            impvec(idx) = -Inf;
        else
            % no singularity : Compute how much improvement was achieved
            if options.verbosity >= 1
                fprintf('\nAccepted split of component %d with cost improvement %g\n', idx, costdiff)
            end
            newNum = newD.num();
            idxOthers = D.invertindex(idx);
            idxNewOthers = idxMap(idxOthers);
            
            % update impvec
            newimpvec = Inf(1, newNum);
            newimpvec(idxNewOthers) = impvec(idxOthers);
            newimpvec(idxSplitted) = costdiff; % this is for both splitted components
            
            % update other variables
            impvec = newimpvec;
            cost = newcost;
            D = newD;
            theta = newtheta;
        end
        if (newNum > options.sm.nummax) 
            break
        end
        if isinf(impvec)
            break
        end
    end

end



function cost = calc_cost(D, theta, data, cost_type, options)
% Calculates the value of the cost
    data = mxe_readdata(data);
    switch cost_type
        case 1
            cost = mxe_costgrad(D, theta, data, options) + ...
                D.BIC(theta, data) / 10 / data.size;
            %cost = - D.ll(theta, data) + D.AICmc(theta, data);
            
        case 2
            cost = mxe_costgrad(D, theta, data, options) + ...
                D.dim() / 10 / data.size;
            % ll may have memory problem!
            %cost = - D.ll(theta, data) + D.dim()/10;
            
    end
end
