%% |mxe_stoppingcriterion|
% *Note:* This is a private function.
%
% Check stopping criteria declared in options
%
% *Syntax*
%
%   stop = mxe_stoppingcriterion(D, theta, options, info, last)
%   stop = mxe_stoppingcriterion(D, theta, options, info, last, first)
%   stop = mxe_stoppingcriterion(D, theta, options, info, last, first, special)
%   [stop, reason] = mxe_stoppingcriterion(...)
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% This is a modified version of the stoppingcriterion function from the
% Manopt toolbox (http://www.manopt.org).
%
% Contributors:
%  Reshad Hosseini
%  Mohamadreza Mash'al
%
% Change log: 
%

function [stop, reason] = mxe_stoppingcriterion(D, theta, options, info, last, first, special)
% |info| is checked no more than the range |info(first)| to |info(last)|   (to prevent checking the stopping criterion that terminated the previous estimation)
% |special| fields contain special objects that may need to stop estimation

    if nargin < 6
        first = 1;
    end
    if nargin < 7
        special = struct;
    end
    
    stop = 0;
    reason = '';
    
    if last < first
        return
    end
    
    stats = info(last);

    % Don't check stopping criteria before the specified minimum iterations
    % has passed
    if isfield(stats, 'iter') && ...
       stats.iter < options.miniter
        return
    end
    
    % Target cost attained
    if isfield(stats, 'cost') && ...
       stats.cost <= options.tolcost
        reason = 'Cost tolerance reached.';
        stop = 1;
        return
    end

    % Target gradient norm attained
    if isfield(stats, 'gradnorm') && ...
       stats.gradnorm < options.tolgradnorm
        reason = 'Gradient norm tolerance reached.';
        stop = 2;
        return
    end

    % Alloted time exceeded
    if isfield(stats, 'time') && ...
       stats.time >= options.maxtime
        reason = 'Max time exceeded.';
        stop = 3;
        return
    end

    % Alloted iteration count exceeded
    if isfield(stats, 'iter') && ...
       stats.iter >= options.maxiter
        reason = 'Max iteration count reached.';
        stop = 4;
        return
    end
    
    % Alloted function evaluation count exceeded
    if isfield(stats, 'costevals') && ...
       stats.costevals >= options.maxcostevals
        reason = 'Maximum number of cost evaluations reached.';
        stop = 5;
        return
    end

    % tolcostdiff
    if isfield(stats, 'cost')
        if last >= first+1 && abs(info(last-1).cost - info(last).cost) <= options.tolcostdiff
            reason = 'Cost difference tolerance reached.';
            stop = 101;
            return
        end
    end
    
    % cross validation
    if options.crossval.enabled
        if isfield(stats, 'iter') && isfield(stats, 'cvCost') && isfield(special, 'cv')
            stopnow = special.cv.stopfun(theta, info, last, first);
            if stopnow
                reason = 'Cross-validation stop triggered.';
                stop = 102;
                return
            end
        end
    end

    % visualization figure closed
    if options.visualization.enabled && options.visualization.stoponclose
        if isfield(special, 'vis')
            stopnow = special.vis.closed();
            if stopnow
                reason = 'Visualization figure closed.';
                stop = 103;
                return
            end
        end
    end

    % Check whether the possibly user defined stopping criterion
    % triggers or not.
    if isfield(options, 'stopfun')
        if nargin(options.stopfun) == 5
            userstop = options.stopfun(D, theta, info, last, first);
        else
            userstop = options.stopfun(D, theta, info, last);
        end
        if userstop
            reason = 'User defined stopfun criterion triggered.';
            stop = 6;
            return
        end
    end
    
end
