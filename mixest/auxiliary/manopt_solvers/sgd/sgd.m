%% |sgd|
% Manifold Stochastic Gradient Descent(SGD) minimization algorithm for Manopt
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Original author: Reshad Hosseini, Aug. 30, 2013.
%
% Change log: 
%   Reshad Hosseini, Jun.26,2013: Improving speed when "transpf" is present 
%

function [x cost info] = sgd(problem, x, options)
% Manifold stochastic gradient descent minimization algorithm for MixEst.
%
% function [x cost info] = sgd(problem)
% function [x cost info] = sgd(problem, x0)
% function [x cost info] = sgd(problem, x0, options)
% function [x cost info] = sgd(problem, [], options)
%
% Apply the SGD minimization algorithm to the problem defined
% in the problem structure, starting at x0 if it is provided (otherwise, at
% a random point on the manifold). To specify options whilst not specifying
% an initial guess, give x0 as [] (the empty matrix).
%
%
% None of the options are mandatory. See the documentation for details.
%
% For input/output descriptions, stopping criteria, help on picking a line
% search algorithm etc, see the help for steepestdescent.
%


    % Verify that the problem description is sufficient for the solver.
    if ~canGetCost(problem)
        warning('manopt:getCost', ...
            'No cost provided. The algorithm will likely abort.');
    end
    if ~canGetGradient(problem)
        warning('manopt:getGradient', ...
            'No gradient provided. The algorithm will likely abort.');
    end

    % Set local defaults here
    localdefaults.sgd.stepsize = 0.001;
    localdefaults.sgd.epoch = 10;

    % Merge global and local defaults, then merge w/ user options, if any.
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);

    % Create a store database
    storedb = struct();

    %timetic = tic();

    % If no initial point x is given by the user, generate one at random.
    if ~exist('x', 'var') || isempty(x)
        x = problem.M.rand();
    end
    
    check_decrease = true;

    if options.verbosity >= 2 || options.sgd.svrg
        if options.verbosity >= 2 && ~ options.sgd.svrg
            warning('Verbosity 2 slows down the procedure');
        end
        [cost egrad storedb] = getEuclideanCostGrad(problem, x, storedb);
        grad= problem.M.egrad2rgrad(x, egrad);
        gradnorm = problem.M.norm(x, grad);
        if options.sgd.svrg
            base_x = x;
            base_egrad = egrad;
            base_grad = grad;
        end
        if options.verbosity >= 2
            if check_decrease
                cost_old = cost;
                x_old = x;
            end
            fprintf(' epoch \t cost val\t grad. norm\n');
            fprintf('%5d\t%+.4e\t%.4e\n', 0, cost, gradnorm);
        end
    end

    epoch = 0; 
    stats = savestats();
    info(1) = stats;
    info(min(10000, options.sgd.epoch+1)).iter = [];
        
    cost = nan;
    % Start iterating until stopping criterion triggers
    for epoch = 1:options.sgd.epoch
        
        % Start timing this iteration
        timetic = tic();
        
        % Random permutation of indicies for each epoch
        data_size = problem.data_size;
        indicies = randperm(data_size);
        batchnum = options.sgd.batchnum;
        batch_size = floor(data_size / batchnum);
        
        for batchIndex = 1:options.sgd.batchnum
            if mod(batchIndex,max(round(options.sgd.batchnum/10),100)) == 1
                fprintf('.')
            end
            
            % selecting some part of indicies
            index_begin = (batchIndex-1) * batch_size + 1;
            index_end = min(batchIndex*batch_size, data_size);
            if batchIndex == batchnum
                index_end = data_size;
            end
            batchIndicies = indicies(index_begin:index_end);
            
            % Compute the new gradient-related quantities for x
            egrad = problem.egradbatch(x, batchIndicies);
            
            
            if length(options.sgd.stepsize) == 1
                if options.sgd.base ~= 1
                    alpha = options.sgd.stepsize * options.sgd.base^((epoch-1)* ...
                            options.sgd.batchnum+batchIndex-1);
                else
                    if ~isinf(options.sgd.diminishc)
                        alpha = options.sgd.stepsize * (options.sgd.diminishc /...
                            ( options.sgd.diminishc + ((epoch-1)* ...
                            options.sgd.batchnum+batchIndex-1).^options.sgd.power ));
                    else
                        alpha = options.sgd.stepsize;
                    end
                end
            else
                lstepsize = length(options.sgd.stepsize);
                nchoose = floor((epoch-1)/options.sgd.epoch*lstepsize)+1;
                alpha = options.sgd.stepsize(nchoose);
            end

            if options.sgd.momentum == 0 || (batchIndex == 1 && epoch == 1)
                % Pick the descent direction as minus the gradient
                desc_dir_euc = problem.D.scaleparam(-1*alpha, egrad);
                desc_dir = problem.M.egrad2rgrad(x, desc_dir_euc);
            else
                if options.sgd.euclidbase
                    % Using Euclidean Gradient with momentum and then retr
                    desc_dir_euc = problem.D.scaleparam(-1*alpha, egrad);
                    desc_dir_euc_old = problem.D.scaleparam(options.sgd.momentum, desc_dir_euc_old);
                    desc_dir_euc = problem.D.sumparam(desc_dir_euc, desc_dir_euc_old);
                    desc_dir = problem.M.egrad2rgrad(x, desc_dir_euc);
                else
                    % Concept of momentum makes sense in Riemmanian domain
                    % Using the same direction for improvement
                    grad = problem.M.egrad2rgrad(x, egrad);
                    distp = problem.M.transp(xold, x, desc_dir_old);
                    desc_dir = problem.M.lincomb(x, -1*alpha, grad, ...
                        options.sgd.momentum, distp);
                end
            end
            if options.sgd.svrg
                if options.sgd.euclidbase
                    % Concept of reducing variance of gradient
                    % We observed it is better to do that on Euclidean domain
                    % Than using Parallel transport in Riemmanian
                    grad_svrg = problem.egradbatch(base_x, batchIndicies);
                    grad_svrg = problem.D.scaleparam(-1, grad_svrg);
                    grad_svrg = problem.D.sumparam(base_egrad, grad_svrg);
                    desc_dir_svrg = problem.D.scaleparam(-1*alpha, grad_svrg);
                    desc_dir_euc = problem.D.sumparam(desc_dir_euc, desc_dir_svrg);
                    desc_dir = problem.M.egrad2rgrad(x, desc_dir_euc);
                else
                    % Different version using parallel transport Riem grad
                    grad_svrg = problem.gradbatch(base_x, batchIndicies);
                    grad_svrg = problem.M.lincomb(base_x, 1, base_grad, -1, grad_svrg);
                    desc_dir_svrg = problem.M.transp(base_x, x, grad_svrg);
                    desc_dir = problem.M.lincomb(x, 1, desc_dir, ...
                        -1*alpha, desc_dir_svrg);
                end
            end
            
            if any(isnan(obj2vec(egrad)))
                break;
            end
            
            xold = x;
            desc_dir_old = desc_dir;
            desc_dir_euc_old = desc_dir_euc;
            x = problem.M.retr(x, desc_dir, 1); 
            
        end
        timetoc = toc(timetic);
        
        % Display iteration information
        if options.verbosity >= 2 || options.sgd.svrg
            % Compute objective-related quantities for x
            [cost egrad storedb] = getEuclideanCostGrad(problem, x, storedb);
            grad = problem.M.egrad2rgrad(x, egrad);
            gradnorm = problem.M.norm(x, grad);
            if options.sgd.svrg
                base_x = x;
                base_egrad = egrad;
                base_grad = grad;
            end
            if options.verbosity >= 2
                if check_decrease
                    if cost > cost_old || any(isnan(obj2vec(egrad)))
                        x = x_old;
                    else
                        cost_old = cost;
                        x_old = x;
                    end
                end
                fprintf('%5d\t%+.4e\t%.4e\n', epoch, cost, gradnorm);
            end
        end
        
        % Save stats in a struct array info, and preallocate
        % (see http://people.csail.mit.edu/jskelly/blog/?x=entry:entry091030-033941)
        stats = savestats();
        info(epoch+1) = stats; %#ok<AGROW>
        
        % Run standard stopping criterion checks
        [stop reason] = stoppingcriterion(problem, x, options, ...
            info, epoch+1);

        if stop
            if options.verbosity >= 1
                fprintf([reason '\n']);
            end
            break;
        end
        
        
    end

    info = info(1:options.sgd.epoch+1);

    if options.verbosity >= 1
        fprintf('Total time is %f [s] (excludes statsfun)\n', ...
                info(end).time);
    end


    % Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = epoch;
        if options.verbosity >= 2
            stats.gradnorm = gradnorm;
            stats.cost = cost;     
        end
        if epoch == 0
            stats.time = 0;%toc(timetic);
        else
            stats.time = info(epoch).time + timetoc;
        end
        
        % Temporary adding
        stats.theta = x;
        stats.ll = -cost;
        
        
        stats = applyStatsfun(problem, x, storedb, options, stats);
    end

end
