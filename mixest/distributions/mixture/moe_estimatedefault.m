%% |moe_estimatedefault|
% *Note:* This is a private function.
%
% Default estimation function for the mixtures of experts. Implements the
% Expectation Maximization (EM) algorithm.
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

function [theta, D, info, options] = moe_estimatedefault(D, data, options)

    data = mxe_readdata(data, false);
    weight = data.weight;
    idxAll = data.index;
    datamat = data.data;

    if nargin < 3
        options = mxe_options();
    else
        options = mxe_options(options);
    end
    
    % cross validation init
    nAll = size(idxAll, 2);
    nTrain = nAll;
    cv = struct;
    if options.crossval.enabled
        [cv, idxTrain] = mxe_crossvalidation(data, ...
            options.crossval.fraction, options.crossval.toliter, ...
            @(theta, data) -llfun(theta, data));
        if ~isempty(weight)
            weight = weight(idxTrain);
        end
        nTrain = cv.nTrain;
    else
        idxTrain = idxAll;
    end

    % plotting init
    if options.plotcost.enabled
        if options.crossval.enabled
            plotCost = mxe_plot({'train-data cost', 'cross-validation cost'}, options.plotcost);
        else
            plotCost = mxe_plot({'cost'}, options.plotcost);
        end
    end
    
    % visualization init
    vis = struct;
    if options.visualization.enabled
        if ~isempty(options.visualization.visobject)
            vis = options.visualization.visobject;
        else
            vis = mxe_visualization(data, options.visualization);
        end
    end
    
    % theta0
    if isempty(options.theta0) 
        if isfield(D, 'init')
            options.theta0 = D.init(struct('data',datamat, 'index',idxTrain));
        else
            options.theta0 = D.randparam();
        end
    end
    theta = options.theta0;

    dataTrain = struct('data',datamat, 'weight',weight, 'index',idxTrain);

    % E-step: calculating the hX parameters
    [llvec, store] = D.llvec(theta, dataTrain);
    hX = exp(bsxfun(@minus, store.logpcond, store.logpcondsum));
    ll = mean(llvec);
    ll_old = ll;   
    ll_diff = 0; 
    
    if isempty(options.previnfo)
        timetic = tic();
        iter = 0;
        first = 1;
        stats = savestats();
        info(1) = stats;
    else
        info = options.previnfo;
        iter = info(end).iter;
        first = numel(info) + 1;
    end
    
    info(min(iter+10000, options.maxiter+1)).iter = [];
    
    while true
        timetic = tic();
        
        % stopping criterion
        [stop, reason] = mxe_stoppingcriterion(D, theta, options, info, iter+1, first, struct('cv',cv, 'vis',vis));
        if stop
            if options.verbosity >= 1
                fprintf([reason '\n']);
            end
            break
        end
            
        
        if options.verbosity >= 2
            fprintf('Log likelihood = %g , LL diff= %g \n', ll, ll_diff);
        end

        % plots update
        plotsUpdated = false;
        if options.plotcost.enabled
            if plotCost.needs_update(iter, 1)
                [ploty, plotx] = mxe_getplotdata('cost', options.plotcost, info, iter+1);
                plotCost.update(plotx, ploty, 1);
                plotsUpdated = true;
            end
            if options.crossval.enabled && plotCost.needs_update(iter, 2)
                [ploty, plotx] = mxe_getplotdata('cvCost', options.plotcost, info, iter+1);
                plotCost.update(plotx, ploty, 2);
                plotsUpdated = true;
            end
        end
        if plotsUpdated
            drawnow
        end
        
        % visualization update
        if options.visualization.enabled
            stop = vis.update(theta, D);
            if stop
                if options.verbosity >= 1
                    fprintf('Visualization figure closed.\n')
                end
                break
            end
        end
        
        % M-step: Update the covariance matrix
        for k = 1:D.num()
            % Inner mixture should run quickly     %TODO ?
%             comp_options = options;
            comp_options = options; 
            comp_options.solver = 'default';
            comp_options.verbosity = 0;
            comp_options.maxiter = 1;
            comp_options.crossval.enabled = false;
            Component = D.expert(k);
            if ~isfield(Component,'estimatedefault');
                comp_options.solver = 'cg';
                comp_options.maxiter = 10;
            end
            % Updating experts
            comp_options.theta0 = theta.D{k};
            theta.D{k} = Component.estimate(struct('data',datamat, 'weight',hX(k,:), 'index',idxTrain), comp_options);           
        end
        % Updating gate
        Gate = D.gate();
        datatot = [datamat(1:Gate.datadim(),idxTrain); hX];
        
        if ~isfield(Component,'estimatedefault');
            comp_options.solver = 'cg';
            comp_options.maxiter = 50;
        end
        comp_options.theta0 = theta.G;
        theta.G = Gate.estimateMstep(datatot, comp_options);
        % E-step: Simultaneous likelihood and weight calculation
        [llvec, store] = D.llvec(theta, dataTrain);
        hX = exp(bsxfun(@minus, store.logpcond, store.logpcondsum));
        ll = mean(llvec);
        ll_diff = ll - ll_old;
        
        ll_old = ll;
        iter = iter + 1;
        stats = savestats();
        info(iter+1) = stats; %#ok<AGROW>
    end
    
    if options.crossval.enabled && cv.hasFinalResult()
        theta = cv.finalResult();
    end
    
    info = info(1:iter+1);

    if options.verbosity >= 1
        fprintf('Total time is %f [s] (excludes statsfun)\n', info(end).time);
    end
    
    
    
    
    function stats = savestats()
        stats.iter = iter;
        stats.cost = -ll;
        if iter == 0
            stats.time = toc(timetic);
        else
            stats.time = info(iter).time + toc(timetic);
        end
        if options.crossval.enabled
            stats = cv.statsfun(theta, stats);
        end
        
        if isfield(options, 'statsfun')
            stats = options.statsfun(D, theta, stats);
        end
    end

    function cost = llfun(theta, data)
        sdata = mxe_readdata(data, false);
        cost = D.ll(theta, data) / sdata.size;
    end
end
