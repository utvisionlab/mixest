%% |smem|
% Split and Merge Expectation-Maximization (SMEM) algorithm for estimation
% of mixture model structure and parameters, with generalized internal
% estimations.
%
% *Syntax*
%
%   [theta, D] = smem(data, num)
%   [theta, D] = smem(data, num, options)
%   [theta, D, info] = smem(...)
%   [theta, D, info, options] = smem(...)
%
% *Description*
%
% This function implements the SMEM algorithm proposed by Ueda et al.
% (2000), though here we use generalized optimization methods.
%
% |[theta, D] = smem(data, num)| returns the estimated parameters and the
% mixture distribution |D| with |num| components, fitted to |data|.
%
% |[theta, D] = smem(data, num, options)| utilizes applicable options from
% the |options| structure in the estimation procedure.
%
% |[theta, D, info] = smem(...)| also returns |info|, a structure array
% containing information about successive iterations performed by iterative
% estimation functions.
%
% |[theta, D, info, options] = smem(...)| also returns the effective
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
% * *|maxCands|* (default 5) : Maximum number of merge candidates.
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
%   options.sm.maxCands = 3;
%   options.inner.partial.maxiter = 10;
%   % fit mixture model to data
%   [theta, D] = smem(data, 5, options)
%
% *References*
%
% # N. Ueda, R. Nakano, Z. Ghahramani, and G. E. Hinton, “Split and Merge
% EM Algorithm for Improving Gaussian Mixture Density Estimates,” The
% Journal of VLSI Signal Processing-Systems for Signal, Image, and Video
% Technology, vol. 26, no. 1–2, pp. 133–140, Aug. 2000.
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

function [theta, D, info, options] = smem(data, num, options)

    assert(num >= 3, 'num should be at least 3.')
    
    extra_options = struct(...
        'sm', struct(...
            'costtype', 2, ...
            'maxcands', 5, ...
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

    % First Run of Full-EM
    if options.verbosity >= 1
        fprintf('\nFirst Run of Full-EM...\n')
    end
    D = mixturefactory(ComponentD, num);
    [theta, D, info] = D.estimate(data, opt_full);
    
    cost_current = calc_cost(D, theta, data, options.sm.costtype);

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


        % Sort the split and merge candidates
        idx_split = D.splitcandidates(theta, data, options);
        [idx_merge1, idx_merge2] = D.mergecandidates(theta, data, options);
        [i, j, k] = combine_candidates(idx_merge1, idx_merge2, idx_split);
        
        maxcands = min(options.sm.maxcands, numel(i));
        
        for c = 1:maxcands

            % check if visualization figure is closed
            if options.visualization.enabled && options.visualization.stoponclose
                if vis.closed()
                    return
                end
            end

            % merge
            if options.verbosity >= 1
                fprintf('\nMerging components (%d, %d) and splitting component (%d)...\n', i(c), j(c), k(c))
            end
            [newD, newtheta, idxMerged, idxMap] = D.merge(i(c), j(c), theta, options, data);
            [newD, newtheta, idxSplitted, idxMap] = newD.split(idxMap(k(c)), newtheta, options, data);
            idxMerged = idxMap(idxMerged);

            % Estimating the parameters partially for the new components
            if options.verbosity >= 1
                fprintf('\nEstimating the parameters partially for the new components...\n')
            end
            opt = mxe_adjustoptions(opt_partial, info(end));
            opt.previnfo = info;
            [newtheta, newD, info] = newD.estimatepartial([idxMerged, idxSplitted(1), idxSplitted(2)], newtheta, data, opt);

            % Full estimation
            if options.verbosity >= 1
                fprintf('\nPerforming full estimation...\n')
            end
            opt = mxe_adjustoptions(opt_full, info(end));
            opt.previnfo = info;
            opt.theta0 = newtheta;
            [newtheta, newD, info] = newD.estimate(data, opt);

            % check improvement
            cost_new = calc_cost(newD, newtheta, data, options.sm.costtype);
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
            cost_current = cost_new;
        else
            if options.verbosity >= 1
                fprintf('\nAll candidates failed to improve the cost. Exitting...\n')
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



function cost = calc_cost(D, theta, data, cost_type)
% Calculates the value of the cost
    switch cost_type
        case 1
            F = D.ll(theta, data) + D.BIC(data);
        case 2
            F = D.ll(theta, data) - D.dim()/100;
    end

    cost = -F;
end



function [i, j, k] = combine_candidates(idx_merge1, idx_merge2, idx_split)
% combine the split and merge candidates (SMEM 2000)

    n_merge = numel(idx_merge1);
    n_split = numel(idx_split);
    
    n_result = n_merge * (n_split-2);
    i = zeros(n_result, 1);
    j = zeros(n_result, 1);
    k = zeros(n_result, 1);
    
    counter = 0;
    for l = 1:n_merge
        il = idx_merge1(l);
        jl = idx_merge2(l);
        for m = 1:n_split
            if idx_split(m)~=il && idx_split(m)~=jl
                counter = counter + 1;
                i(counter) = il;
                j(counter) = jl;
                k(counter) = idx_split(m);
            end
        end
    end
end
