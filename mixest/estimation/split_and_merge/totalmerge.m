%% |totalmerge|
% Total-merge algorithm for estimation of mixture model structure and
% parameters, with generalized internal estimations. This function is run
% after running a function like smem or totalsplit
%
% *Syntax*
%
%   [theta, D, info, options] = totalmerge(data, theta0, target_num, options)
%
% *Description*
%
% This function implements the total-merge algorithm -- It is for merging
% candidates based on certain criteria to get the final results for 
% different number of components
%
% |totalmerge(data, theta0, target_num,  options)| estimates and saves 
% parameters and the mixture distributions |D| from number of
% components in |theta0| to |target_num| components,fitted to |data|. It
% can also return the last merge stage
% Point: If info for all merges is needed one should do only one merge
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
% * *|savefoler|* (default result) : folder to save the results for
% different components -- the files are called tmerge%i where %i is number 
% of components
%
% |options.inner| can be used to set specific options for the inner
% estimations. To set options only for the partial or full inner
% estimations use |options.inner.partial| or |options.inner.full|
% respectively.
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

function [theta, D, info, options] = totalmerge(data, theta0, target_num, options)
% This function implements the Total-Merge method 
%
% This procedure applies merging approach (partial merging to increase
% the speed). 
%
    newNum = length(theta0.D);
    
    extra_options = struct(...
        'sm', struct(...
            'costtype', 1, ...
            'tolcostdiff', 0, ...
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
        D = mixturefactory(ComponentD, length(theta0.D));
    else
        ds = mxe_readdata(data, false);
        ComponentD = mvnfactory(ds.dim);
        D = mixturefactory(ComponentD, length(theta0.D));
        clear ds
    end
    
    % Creating a directory for saving the results
    mkdir(options.sm.savefolder);
    
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
        fprintf('\nStarting with newNum component...\n')
    end
    opt_full.theta0 = theta0;
    [theta, D, info] = D.estimate(data, opt_full);
    
    namesave = [options.sm.savefolder '/tmerge' num2str(newNum)];
    save(namesave, 'theta', 'D');

    % Merging loop
    while D.num() > target_num
        
        % check if visualization figure is closed
        if options.visualization.enabled && options.visualization.stoponclose
            if vis.closed()
                return
            end
        end
        
        newNum = newNum - 1;
        
        % merge the best candidates
        [idx1, idx2] = D.mergecandidates(theta, data, options, 1);
        if options.verbosity >= 1
            fprintf('\nMerging components (%d, %d) from the mixture with %d components...\n', idx1, idx2, D.num())
        end
        [D, theta, idxMerged] = D.merge(idx1, idx2, theta, options, data);
        
        % Estimating the parameters partially for the newly merged component
        [theta, D] = D.estimatepartial(idxMerged, theta, data, opt_partial);
        
        opt_full.theta0 = theta;
        [theta, D, info] = D.estimate(data, opt_full);
        
        namesave = [options.sm.savefolder '/tmerge' num2str(newNum)];
        save(namesave, 'theta', 'D');
        
    end

end