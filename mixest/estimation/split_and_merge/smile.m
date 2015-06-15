%% |smile|
% Split–Merge Incremental LEarning (SMILE) algorithm for estimation
% of mixture model structure and parameters, with generalized internal
% estimations.
%
% *Syntax*
%
%   [theta, D] = smile(data, target_num)
%   [theta, D] = smile(data, target_num, options)
%   [theta, D, info] = smile(...)
%   [theta, D, info, options] = smile(...)
%
% *Description*
%
% This function implements the SMILE algorithm proposed by Blekas and
% Lagaris (2007), though here we use generalized optimization methods.
%
% |[theta, D] = smile(data, target_num)| returns the estimated parameters
% and the mixture distribution |D| with |target_num| components, fitted to
% |data|.
%
% |[theta, D] = smile(data, target_num, options)| utilizes applicable
% options from the |options| structure in the estimation procedure.
%
% |[theta, D, info] = smile(...)| also returns |info|, a structure array
% containing information about successive iterations performed by iterative
% estimation functions.
%
% |[theta, D, info, options] = smile(...)| also returns the effective
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
% * *|numInit|* (default |2|) : Initial number of mixture components.
% * *|tolCostDiff|* (default |0|) : Minimum decrease in cost to accept a
% SOMO operation.
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
% |info| returned by the Manopt solver is returned directly. For the
% 'default' solver see the documentation of the 'estimatedefault' function
% for the mixture distribution. You can read more at our documentation on
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
%   options.sm.numInit = 3;
%   options.inner.partial.maxiter = 10;
%   % fit mixture model to data
%   [theta, D] = smile(data, 5, options)
%
% *References*
%
% # K. Blekas and I. E. Lagaris, “Split–Merge Incremental LEarning (SMILE)
% of Mixture Models,” in Artificial Neural Networks – ICANN 2007, J. M. de
% S?, L. A. Alexandre, W. Duch, and D. Mandic, Eds. Springer Berlin
% Heidelberg, 2007, pp. 291–300.

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Contributors:
%  Reshad Hosseini
%  Mohamadreza Mash'al
%
% Change log: 
%

function [theta, D, info, options] = smile(data, target_num, options)

    extra_options = struct(...
        'sm', struct(...
            'numinit', 2, ...
            'tolcostdiff', 0, ...
            'componentd', [] ... default is mvnfactory(datadim)
        ));

    if nargin < 2
        options = mxe_options([], extra_options);
    else
        options = mxe_options(options, extra_options);
    end
        
    if isempty(options.sm.componentd)
        ds = mxe_readdata(data, false);
        ComponentD = mvnfactory(ds.dim);
        clear ds
    else
        ComponentD = options.sm.componentd;
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

    % construct a mixture with two components
    k = options.sm.numinit;
    if options.verbosity >= 1
        fprintf('\nStarting with %d components...\n', k)
    end
    D = mixturefactory(ComponentD, k);
    [theta, D, info] = D.estimate(data, opt_full);
    
    % Estimate the current cost
    cost = calc_cost(D, theta, data);
    
    while k < target_num
        
        % check if visualization figure is closed
        if options.visualization.enabled && options.visualization.stoponclose
            if vis.closed()
                return
            end
        end
        
        % Split: select a cluster and divide it into two clusters
        idx = D.splitcandidates(theta, data, options, 1);
        if options.verbosity >= 1
            fprintf('\nSplitting component %d from %d mixture components...\n', idx, D.num())
        end
        [newD, newtheta, idxSplitted] = D.split(idx, theta, options, data);
        
        % Optimization operation: Perform partial-EM and then full EM
        if options.verbosity >= 1
            fprintf('\nPerforming partial estimation...\n')
        end
        opt = mxe_adjustoptions(opt_partial, info(end));
        opt.previnfo = info;
        [newtheta, newD, info] = newD.estimatepartial(idxSplitted, newtheta, data, opt);
        if options.verbosity >= 1
            fprintf('\nPerforming full estimation...\n')
        end
        opt = mxe_adjustoptions(opt_full, info(end));
        opt.previnfo = info;
        opt.theta0 = newtheta;
        [newtheta, newD, info] = newD.estimate(data, opt);
        
        % split operation is always accepted
        theta = newtheta;
        D = newD;
        
        
        % Merge: select two clusters and merge them
        [idx1, idx2] = D.mergecandidates(theta, data, options, 1);
        if options.verbosity >= 1
            fprintf('\nMerging components (%d, %d) from the mixture with %d components...\n', idx1, idx2, D.num())
        end
        [newD, newtheta, idxMerged] = D.merge(idx1, idx2, theta, options, data);
        
        % Optimization operation: Perform partial-EM and then full EM
        if options.verbosity >= 1
            fprintf('\nPerforming partial estimation...\n')
        end
        opt = mxe_adjustoptions(opt_partial, info(end));
        opt.previnfo = info;
        [newtheta, newD, info] = newD.estimatepartial(idxMerged, newtheta, data, opt);
        if options.verbosity >= 1
            fprintf('\nPerforming full estimation...\n')
        end
        opt = mxe_adjustoptions(opt_full, info(end));
        opt.previnfo = info;
        opt.theta0 = newtheta;
        [newtheta, newD, info] = newD.estimate(data, opt);
        
        
        
        % calculate the cost improvement
        newcost = calc_cost(newD, newtheta, data);
        costdiff = cost - newcost;
        
        if costdiff < options.sm.tolcostdiff
            if options.verbosity >= 1
                fprintf('\nRejected the last merge operation. Cost improvement %g\n', costdiff)
            end
            k = k + 1;
        else
            if options.verbosity >= 1
                fprintf('\nAccepted the fused cluster. Cost improvement %g\n', costdiff)
            end
            cost = newcost;
            theta = newtheta;
            D = newD;
        end
    end
    
    if options.verbosity >= 1
        fprintf('\nReached the target number of components (%d). Exitting...\n', target_num)
    end
    
    % visualization update
    if options.visualization.enabled
        vis.update(theta, D);
    end
end


function cost = calc_cost(D, theta, data)
% Calculates the value of the cost

    F = D.ll(theta, data) + D.BIC(data);
    cost = -F;
end

