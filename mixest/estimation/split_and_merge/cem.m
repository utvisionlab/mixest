%% |cem|
% Competitive Expectation-Maximization (CEM) algorithm for estimation
% of mixture model structure and parameters, with generalized internal
% estimations.
%
% *Syntax*
%
%   [theta, D] = cem(data)
%   [theta, D] = cem(data, options)
%   [theta, D, info] = cem(...)
%   [theta, D, info, options] = cem(...)
%
% *Description*
%
% This function implements the CEM algorithm proposed by Zhang et al.
% (2004), though here we use generalized optimization methods.
%
% |[theta, D] = cem(data)| returns the estimated parameters and the
% mixture distribution |D|, fitted to |data|.
%
% |[theta, D] = cem(data, options)| utilizes applicable options from the
% |options| structure in the estimation procedure.
%
% |[theta, D, info] = cem(...)| also returns |info|, a structure array
% containing information about successive iterations performed by iterative
% estimation functions.
%
% |[theta, D, info, options] = cem(...)| also returns the effective
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
% * *|numInit|* (default 1) : Initial number of mixture components.
% * *|numMin|* (default 1) : Minimum number of mixture components.
% * *|numMax|* (default 15) : Maximum number of mixture components.
% * *|mergeMaxCands|* (default 5) : Maximum number of merge candidates.
% * *|splitMaxCands|* (default 5) : Maximum number of split candidates.
% * *|maxFail|* (default 2) : Maximum number of split/merge trials before exit.
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
%   options.sm.numMax = 7;
%   options.inner.partial.maxiter = 10;
%   % fit mixture model to data
%   [theta, D] = cem(data, options)
%
% *References*
%
% # B. Zhang, C. Zhang, and X. Yi, “Competitive EM algorithm for finite
% mixture models,” Pattern Recognition, vol. 37, no. 1, pp. 131–144, Jan.
% 2004.
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% This code is inspired by the FSMEM_MVGM code by Sebastien Paris.
%
% Contributors:
%  Reshad Hosseini
%  Mohamadreza Mash'al
%
% Change log: 
%

function [theta, D, info, options] = cem(data, options)

    extra_options = struct(...
        'sm', struct(...
            'numinit', 1, ...
            'nummin', 1, ...
            'nummax', 15, ...
            'mergemaxcands', 5, ...
            'splitmaxcands', 5, ...
            'maxfail', 2, ...
            'componentd', [] ... default is mvnfactory(datadim)
        ));

    if nargin < 2
        options = mxe_options([], extra_options);
    else
        options = mxe_options(options, extra_options);
    end
    
    
    Kini = options.sm.numinit;
    Kmin = options.sm.nummin;
    Kmax = options.sm.nummax;
    maxcands_merge = options.sm.mergemaxcands;
    maxcands_split = options.sm.splitmaxcands;
    fail_exit = options.sm.maxfail;
    
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

    % init
    Kcurrent = Kini;
    fail = 0;
    next = -1; % merge
    
    % First Run of Full-EM
    if options.verbosity >= 1
        fprintf('\nFirst Run of Full-EM. K=%d...\n', Kcurrent)
    end
    D = mixturefactory(ComponentD, Kcurrent);
    [theta, D, info] = D.estimate(data, opt_full);
    
    cost_current = calc_cost(D, theta, data);

    if fail_exit > 0
        if Kcurrent == 1
            next = 1;
        end

        iter = 0;
        while true

            % check special stopping criteria
            if iter >= options.miniter
                
                % check maxiter for the main loop
                if iter >= options.maxiter
                    if options.verbosity >= 1
                        fprintf('Max iteration count for the main loop reached.\n');
                    end
                    break
                end
                
            end
            
            
            if (next < 0) && (Kcurrent > Kmin - 1) % Merge Move
                Knew = Kcurrent - 1;
                cost_new = inf;
                
                [idx_merge1, idx_merge2] = D.mergecandidates(theta, data, options, maxcands_merge);
                
                for m = 1:numel(idx_merge1)
                    % check if visualization figure is closed
                    if options.visualization.enabled && options.visualization.stoponclose
                        if vis.closed()
                            return
                        end
                    end
            
                    idx1 = idx_merge1(m);
                    idx2 = idx_merge2(m);
                    
                    % merge
                    if options.verbosity >= 1
                        fprintf('\nMerging components (%d, %d) from the mixture with %d components...\n', idx1, idx2, D.num())
                    end
                    [newD, newtheta, idxMerged] = D.merge(idx1, idx2, theta, options, data);

                    % Estimating the parameters partially for the newly merged component
                    if options.verbosity >= 1
                        fprintf('\nEstimating the parameters partially for the newly merged component...\n')
                    end
                    opt = mxe_adjustoptions(opt_partial, info(end));
                    opt.previnfo = info;
                    [newtheta, newD, info] = newD.estimatepartial(idxMerged, newtheta, data, opt);
                    
                    % Full estimation
                    if options.verbosity >= 1
                        fprintf('\nPerforming full estimation...\n')
                    end
                    opt = mxe_adjustoptions(opt_full, info(end));
                    opt.previnfo = info;
                    opt.theta0 = newtheta;
                    [newtheta, newD, info] = newD.estimate(data, opt);
                    
                    % check improvement
                    cost_new = calc_cost(newD, newtheta, data);
                    if cost_new < cost_current
                        if options.verbosity >= 1
                            fprintf('\nCost improved :)\n')
                        end
                        break
                    end
                    if options.verbosity >= 1
                        fprintf('\nNo cost improvement :(\n')
                    end
                end
                if cost_new < cost_current
                    theta = newtheta;
                    D = newD;
                    Kcurrent = Knew;
                    cost_current = cost_new;
                    fail = 0;
                    if Kcurrent == Kmin
                        next = 1;
                    end
                else
                    if options.verbosity >= 1
                        fprintf('All merge candidates failed to improve the cost.\n')
                    end
                    fail = fail + 1;
                    next = 1;
                end
                
            elseif (next > 0) && (Kcurrent < Kmax) % Split Move
            
                Knew = Kcurrent + 1;
                cost_new = inf;
                
                idx_split = D.splitcandidates(theta, data, options, maxcands_split);
                
                for m = 1:numel(idx_split)
                    % check if visualization figure is closed
                    if options.visualization.enabled && options.visualization.stoponclose
                        if vis.closed()
                            return
                        end
                    end
            
                    idx = idx_split(m);
                    
                    % based on the index split the component
                    if options.verbosity >= 1
                        fprintf('\nSplitting component %d from %d mixture components...\n', idx, D.num())
                    end
                    [newD, newtheta, idxSplitted] = D.split(idx, theta, options, data);

                    % Estimating the parameters partially for the newly splitted components
                    if options.verbosity >= 1
                        fprintf('\nEstimating the parameters partially for the newly splitted components...\n')
                    end
                    opt = mxe_adjustoptions(opt_partial, info(end));
                    opt.previnfo = info;
                    [newtheta, newD, info] = newD.estimatepartial(idxSplitted, newtheta, data, opt);
                    
                    % Full estimation
                    if options.verbosity >= 1
                        fprintf('\nPerforming full estimation...\n')
                    end
                    opt = mxe_adjustoptions(opt_full, info(end));
                    opt.previnfo = info;
                    opt.theta0 = newtheta;
                    [newtheta, newD, info] = newD.estimate(data, opt);
                    
                    % check improvement
                    cost_new = calc_cost(newD, newtheta, data);
                    if cost_new < cost_current
                        if options.verbosity >= 1
                            fprintf('\nCost improved :)\n')
                        end
                        break
                    end
                    if options.verbosity >= 1
                        fprintf('\nNo cost improvement :(\n')
                    end
                end
                if cost_new < cost_current
                    theta = newtheta;
                    D = newD;
                    Kcurrent = Knew;
                    cost_current = cost_new;
                    fail = 0;
                    if Kcurrent == Kmin
                        next = -1;
                    end
                else
                    if options.verbosity >= 1
                        fprintf('All split candidates failed to improve the cost.\n')
                    end
                    fail = fail + 1;
                    next = -1;
                end
                
            elseif (fail > 0) && (Kcurrent == Kmin) || (Kcurrent == Kmax)
                fail = fail_exit;
            end
            
            if fail >= fail_exit
                if options.verbosity >= 1
                    fprintf('\nMaximum failure count reached. Exitting...\n')
                end
                
                % visualization update
                if options.visualization.enabled
                    vis.update(theta, D);
                end
                
                break
            end
            
            
            iter = iter + 1;
        end
    end

end



function cost = calc_cost(D, theta, data)
% Calculates the value of the cost

    F = D.ll(theta, data) + D.BIC(data);
    cost = -F;
end
