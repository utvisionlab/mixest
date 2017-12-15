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

    if options.verbosity >= 2 || options.sgd.svrg
        if options.verbosity >= 2 && ~ options.sgd.svrg
            warning('Verbosity 2 slows down the procedure');
        end
        [cost grad storedb] = getCostGrad(problem, x, storedb);
        gradnorm = problem.M.norm(x, grad);
        if options.sgd.svrg
            base_x = x;
            base_grad = grad;
        end
        if options.verbosity >= 2
            fprintf(' epoch \t cost val\t grad. norm\n');
            fprintf('%5d\t%+.4e\t%.4e\n', 0, cost, gradnorm);
        end
    end

    epoch = 0; 
    stats = savestats();
    info(1) = stats;
    info(min(10000, options.sgd.epoch+1)).iter = [];
        
    % Start iterating until stopping criterion triggers
    for epoch = 1:options.sgd.epoch
        
        % Start timing this iteration
        timetic = tic();

        for batchIndex = 1:options.sgd.batchnum
            
            % Compute the new gradient-related quantities for x
            grad = problem.gradbatch(x, batchIndex);
            
            
            if ~isinf(options.sgd.diminishC)
                alpha = options.sgd.stepsize * (options.sgd.diminishC /...
                    (options.sgd.diminishC + (epoch-1)*options.sgd.batchnum+batchIndex-1));
            else
                alpha = options.sgd.stepsize;
            end
            if options.sgd.momentum == 0 || (batchIndex == 1 && epoch == 1)
                % Pick the descent direction as minus the gradient
                desc_dir = problem.M.lincomb(x, -1*alpha, grad);
            else
                distp = problem.M.log(x, xold);
                desc_dir = problem.M.lincomb(x, -1*alpha, grad, ...
                    -1*options.sgd.momentum, distp);
            end
            if options.sgd.svrg
                grad_svrg = problem.gradbatch(base_x, batchIndex);
                grad_svrg = problem.M.lincomb(base_x, 1, base_grad, -1, grad_svrg);
                desc_dir_svrg = problem.M.transp(base_x, x, grad_svrg);
                desc_dir = problem.M.lincomb(x, 1, desc_dir, ...
                    -1*alpha, desc_dir_svrg);
            end
            xold = x;
            x = problem.M.retr(x, desc_dir, 1); 
            
        end
        
        timetoc = toc(timetic);
        
        % Display iteration information
        if options.verbosity >= 2 || options.sgd.svrg
            % Compute objective-related quantities for x
            [cost grad storedb] = getCostGrad(problem, x, storedb);
            gradnorm = problem.M.norm(x, grad);
            if options.sgd.svrg
                base_x = x;
                base_grad = grad;
            end
            if options.verbosity >= 2
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
        stats = applyStatsfun(problem, x, storedb, options, stats);
    end

end
