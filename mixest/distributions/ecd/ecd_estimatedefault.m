%% |ecd_estimatedefault|
% *Note:* This is a private function.
%
% Default estimation function for ecd distribution. Implements alternative
% coordinate descent (between Radial and Sigma).
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

function [theta, D, info, options] = ecd_estimatedefault(D, data, options)

    data = mxe_readdata(data, false);
    weight = data.weight;
    idxAll = data.index;
    datamat = data.data;

    if nargin < 3
        options = mxe_options();
    else
        options = mxe_options(options);
    end
    
    % penalizer_theta
    if options.penalize
        if isempty(options.penalizertheta)
            options.penalizertheta = D.penalizerparam(data);
        end
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

    dataRadial = struct('data',zeros(1,size(datamat,2)), ...
        'weight',weight, 'index',idxTrain);
    radialD = D.radialD();
    
    % Computing Square Radial and ll
    [u, logdetR, radialDllvec]  = compute_radial(theta, datamat);
    llvec = - (0.5*datadim)*log(pi) - logdetR + gammaln(0.5*datadim) ...
        + radialDllvec + (1-0.5*datadim)*log(u);
    if ~isempty(weight)
        ll = mean(llvec(idxTrain) .* weight(idxTrain), 2);
    else
        ll = mean(llvec(idxTrain), 2);
    end
    %ll
    
    %sum(D.llvec(theta, dataTrain))/nTrain
    if options.penalize
        costPen = D.penalizercost(theta, options.penalizertheta);
        ll = ll + costPen / length(llvec);
    end
    
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
        [stop, reason] = mxe_stoppingcriterion(D, theta, options, info, ...
            iter+1, first, struct('cv',cv, 'vis',vis));
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
        
        % Update the parameters
        comp_options = options;
        comp_options.solver = 'default';
        comp_options.verbosity = 0;
        comp_options.miniter = 2;
        comp_options.maxiter = 1;
        comp_options.crossval.enabled = false;
        if ~isfield(radialD,'estimatedefault');
            comp_options.solver = 'cg';
            comp_options.maxiter = 10;
        end
        if options.penalize
            comp_options.penalizertheta = options.penalizertheta.radialD;
        end
        comp_options.theta0 = theta.radialD;
        dataRadial.data = u;
        thetaRadial = radialD.estimatedefault(dataRadial, comp_options);
        comp_options.maxiter = 1;
        comp_options.theta0 = theta;
        if strcmp (radialD.name(), 'gamma')
            theta = ecd_eg_estimatedefault(D, thetaRadial, dataTrain, comp_options);
        end
        theta.radialD = thetaRadial;
        if iter < options.maxiter - 1
            % Computing Square Radial and ll
            [u, logdetR, radialDllvec]  = compute_radial(theta, datamat);
            llvec = - (0.5*datadim)*log(pi) - logdetR + gammaln(0.5*datadim) ...
                + radialDllvec + (1-0.5*datadim)*log(u);
            if ~isempty(weight)
                ll = mean(llvec(idxTrain) .* weight(idxTrain), 2);
            else
                ll = mean(llvec(idxTrain), 2);
            end
            if options.penalize
                costPen = D.penalizercost(theta, options.penalizertheta);
                ll = ll + costPen / length(llvec);
            end
        end
        
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

    function [u, logdetR, radialDllvec]  = compute_radial(theta, data)
        data = mxe_readdata(data);
        data = data.data;
        datadim = size(data,1);
        [R, p] = chol(theta.sigma); % C=R' R ( R upper trangular)
        if p > 0
            warning('matrix is not symmetric positive definite ...');
            t = 1e-10;
            while p > 0
                theta.sigma = theta.sigma + t * eye(datadim);
                [R, p] = chol(theta.sigma); %#ok<ASGLU>
                t = t * 100;
            end
            R = chol(theta.sigma);
        end
        logdetR = sum(log(diag(R)));
        Rinv = R \ eye(datadim); % Faster version of Rinv = inv_triu(R);
        % u = X' C^-1 X
        Rinvdata = (Rinv' * data);
        u = sum(Rinvdata.^2, 1);
        radialDllvec = radialD.llvec(theta.radialD, u);
    end
end
            
