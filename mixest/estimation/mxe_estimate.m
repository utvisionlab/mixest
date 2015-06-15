%% |mxe_estimate|
% *Note:* This is a private function.
%
% Main estimation function using manifold optimization
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

function [theta, D, info, options] = mxe_estimate(D, data, options)

    if nargin < 3
        options = mxe_options();
    else
        options = mxe_options(options);
    end

    % check if default estimation method should be used
    if strcmpi(options.solver, 'default') && isfield(D, 'estimatedefault')
        if options.verbosity >= 1
            fprintf('Using default estimation algorithm...\n')
        end
        [theta, D, info, options] = D.estimatedefault(data, options);
        return
    end

    data = mxe_readdata(data, false);
    idxAll = data.index;
    datamat = data.data;

    % get solver handle (if it is not already)
    if ischar(options.solver)
        if strcmpi(options.solver, 'default')
            if options.verbosity >= 1
                fprintf('Using LBFGS solver...\n')
            end
            options.solver = 'lbfgs';
        end
        solver = mxe_getsolverhandle(options.solver);
    else
        solver = options.solver;
    end

    % manipulate required solver options
    if isfield(options, 'stopfun')
        user_stopfun = options.stopfun;
    else
        user_stopfun = [];
    end
    options.stopfun = @stopfun;
    if isfield(options, 'statsfun')
        user_statsfun = options.statsfun;
    else
        user_statsfun = [];
    end
    options.statsfun = @statsfun;


    % penalizer_theta
    if options.penalize
        if isempty(options.penalizertheta)
            options.penalizertheta = D.penalizerparam(data);
        end
    end

    % cross validation init
    if options.crossval.enabled
        [cv, idxTrain] = mxe_crossvalidation(struct('data',datamat, 'index',idxAll), ...
            options.crossval.fraction, options.crossval.toliter, ...
            @(theta, data) mxe_costgrad(D, theta, data, options));
    else
        idxTrain = idxAll;
    end

    % mini batch init
    if options.minibatch.size > 0
        [mb, idxTrain] = mxe_minibatch(idxTrain, options);
    end

    % problem
    %dataTrain = struct('data',datamat, 'index',idxTrain);
    problem.M = D.M;
    if isfield(options, 'costgrad')
        problem.costgrad = options.costgrad;
    else
        problem.costgrad = @costgrad;
    end
    function [cost, grad] = costgrad(theta)
        if nargout > 1
            [cost, grad] = mxe_costgrad(D, theta, ...
                struct('data',datamat, 'index',idxTrain), options);
        else
             cost = mxe_costgrad(D, theta, ...
                 struct('data',datamat, 'index',idxTrain), options);
        end
    end
    % Adding store in costgrad decreases the speed in line search
    %problem.costgrad = @(theta, store) mxe_costgrad(D, theta, dataTrain, options, store);

    
    % checkgradient
    if options.checkgradient
        checkgradient(problem);
        return
    end

    % plotting init
    if options.plotcost.enabled
        if options.crossval.enabled
            plotCost = mxe_plot({'train-data cost', 'cross-validation cost'}, options.plotcost);
        else
            plotCost = mxe_plot({'cost'}, options.plotcost);
        end
    end
    if options.plotgradnorm.enabled
        reset_axes = (options.plotgradnorm.axes ~= options.plotcost.axes);
        plotGradNorm = mxe_plot({'gradient norm'}, options.plotgradnorm, reset_axes);
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

    % theta0
    if isempty(options.theta0) && isfield(D, 'init')
        options.theta0 = D.init(struct('data',datamat, 'index',idxTrain));
    end % else Manopt will generate a random theta0 if it is empty
    
    % TODO: Update manopt optimization algorithms of takig curiter into account
    % iterations passed during current call (Note: We set info(last).iter to
    % the complete iterations number in statsfun)
    curiter = 0;
    % initialize some variables used when calling the solver multiple times
    % to discard optimization history.
    prevtheta = options.theta0; % last theta from the previous call, to be used as the new initial point
    if isempty(options.previnfo)
        allinfo = [];
        allinfolast = 0;
        previter = 0; % previous iterations
        prevtime = 0; % previous execution time
    else
        allinfo = options.previnfo;
        allinfolast = numel(allinfo);
        previter = allinfo(allinfolast).iter; % previous iterations
        prevtime = allinfo(allinfolast).time; % previous execution time
    end
    mbContinue = true; % flag to indicate that we should call the solver again

    % solve
    while mbContinue
        mbContinue = false;
        [theta, cost, callinfo] = solver(problem, prevtheta, options); %#ok<ASGLU>

        [allinfo, allinfolast] = mxe_mergeinfo(allinfo, allinfolast, callinfo, options.maxiter);
    end

    if options.crossval.enabled && cv.hasFinalResult()
        theta = cv.finalResult();
    end
    info = allinfo(1:allinfolast);

    if nargout > 3
        % restore manipulated options
        if isempty(user_stopfun)
            options = rmfield(options, 'stopfun');
        else
            options.stopfun = user_stopfun;
        end
        if isempty(user_statsfun)
            options = rmfield(options, 'statsfun');
        else
            options.statsfun = user_statsfun;
        end
    end

    return






%% nested functions

    function stats = statsfun(problem, theta, stats) %#ok<INUSL>
        
        if previter > 0
            if curiter == 0
                stats.time = stats.time + prevtime;
            end
            stats.iter = stats.iter + previter;
        end
        
        curiter = curiter + 1;
        
        if options.minibatch.size > 0
            if curiter > options.minibatch.iter
                curiter = 0;
                idxTrain = mb.next();
            end      
            
            stats.curiterratio = curiter / options.minibatch.iter;
        end
        
        % If cross-validation is used then also store these values in info
        if options.crossval.enabled
            stats = cv.statsfun(theta, stats);
        end
                
        % if user has specified a statsfun, call it
        if ~isempty(user_statsfun)
            stats = user_statsfun(D, theta, stats);
        end
    end



    function stopnow = stopfun(problem, theta, info, last) %#ok<INUSL>
        
        stopnow = false;
        
        iter = info(last).iter; % Note: iter starts from zero
        % iter: total iterations; curiter: current run iterations  (see statsfun)
        
        % check stopping criteria that Manopt doesn't check
        if iter >= options.miniter
        
            % cross validation
            if options.crossval.enabled
                stopnow = cv.stopfun(theta, info, last);
                if stopnow
                    if options.verbosity >= 1
                        fprintf('Cross-validation stop triggered.\n')
                    end
                    return
                end
            end

            % tolcostdiff
            if isfield(info, 'cost')
                if last >= 2 && abs(info(last-1).cost - info(last).cost) <= options.tolcostdiff
                    if options.verbosity >= 1
                        fprintf('Cost difference tolerance reached.\n')
                    end
                    stopnow = true;
                    return
                end
            end
            
        end % stopping criteria
        
        
        % mini batch
        if options.minibatch.size > 0
            % first iteration of new minibatch is the same as the
            % last iteration of old minibatch.
            if curiter == 0 && options.minibatch.discardhistory
                if options.verbosity >= 2
                    fprintf('%d iterations passed. Moving to next mini batch...\n', iter)
                end
                stopnow = true;
                mbContinue = true;
                previter = iter;
                prevtime = info(last).time;
                prevtheta = theta;
                return
                
            end
        end
        
        % plots update
        plotsUpdated = false;
        if options.plotcost.enabled
            if plotCost.needs_update(iter, 1)
                [ploty, plotx] = mxe_getplotdata('cost', options.plotcost, info, last, previter, allinfo, allinfolast);
                plotCost.update(plotx, ploty, 1);
                plotsUpdated = true;
            end
            if options.crossval.enabled && plotCost.needs_update(iter, 2)
                [ploty, plotx] = mxe_getplotdata('cvCost', options.plotcost, info, last, previter, allinfo, allinfolast);
                plotCost.update(plotx, ploty, 2);
                plotsUpdated = true;
            end
        end
        if options.plotgradnorm.enabled
            if plotGradNorm.needs_update(iter)
                [ploty, plotx] = mxe_getplotdata('gradnorm', options.plotgradnorm, info, last, previter, allinfo, allinfolast);
                plotGradNorm.update(plotx, ploty);
                plotsUpdated = true;
            end
        end
        if plotsUpdated
            drawnow
        end
        
        % visualization update
        if options.visualization.enabled
            stopnow = vis.update(theta, D);
            if stopnow && options.verbosity >= 1
                fprintf('Visualization figure closed.\n')
            end
        end
        
        % if user has specified a stopfun, call it
        if ~stopnow && ~isempty(user_stopfun)
            if previter > 0
                stopnow = user_stopfun(D, theta, [allinfo(1:allinfolast) info(2:last)], allinfolast+last-1);
            else
                stopnow = user_stopfun(D, theta, info, last);
            end
        end
    end

end
