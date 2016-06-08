%% |lbfgs|
% Manifold LBFGS minimization algorithm for Manopt
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Original author: Reshad Hosseini, Aug. 30, 2013.
%
% Change log: 
%   Reshad Hosseini, Jun.26,2013: Improving speed when "transpf" is present 
%

function [x cost info] = lbfgs(problem, x, options)
% Manifold LBFGS minimization algorithm for Manopt.
%
% function [x cost info] = lbfgs(problem)
% function [x cost info] = lbfgs(problem, x0)
% function [x cost info] = lbfgs(problem, x0, options)
% function [x cost info] = lbfgs(problem, [], options)
%
% Apply the LBFGS minimization algorithm to the problem defined
% in the problem structure, starting at x0 if it is provided (otherwise, at
% a random point on the manifold). To specify options whilst not specifying
% an initial guess, give x0 as [] (the empty matrix).
%
% For more information about the algorithm see:
%  Qi et al; ???? 2010. Riemannian BFGS Algorithm with Applications 
%
%
% None of the options are mandatory. See the documentation for details.
%
% For input/output descriptions, stopping criteria, help on picking a line
% search algorithm etc, see the help for steepestdescent.
%
% See also: steepestdescent linesearch
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
localdefaults.minstepsize = 1e-10;
localdefaults.maxiter = 1000;
localdefaults.tolgradnorm = 1e-6;
localdefaults.numgrad = 20;
localdefaults.linesearch = @linesearch_wolfe; %_adaptive
costevals = 1;

% Merge global and local defaults, then merge w/ user options, if any.
localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
if ~exist('options', 'var') || isempty(options)
    options = struct();
end
options = mergeOptions(localdefaults, options);

% Create a store database
storedb = struct();

timetic = tic();

% If no initial point x is given by the user, generate one at random.
if ~exist('x', 'var') || isempty(x)
    x = problem.M.rand();
end

% Compute objective-related quantities for x
[cost grad storedb] = getCostGrad(problem, x, storedb);
gradnorm = problem.M.norm(x, grad);

% Iteration counter (at any point, iter is the number of fully executed
% iterations so far)
iter = 0;

% Save stats in a struct array info, and preallocate
% (see http://people.csail.mit.edu/jskelly/blog/?x=entry:entry091030-033941)
stats = savestats();
info(1) = stats;
info(min(10000, options.maxiter+1)).iter = [];

% Initial line search memory
lsmem = [];

if options.verbosity >= 2
    fprintf(' iter\t func eval \t cost val\t grad. norm\n');
end

% Compute a normalized descent direction
desc_dir = problem.M.lincomb(x, -1/gradnorm, grad);

% Initialize the Hessian
H = 1;
x_all = {x};
desc_dir_all = {};
grad_diff_all = {};
ddgd_all = {};
gd_all = {};
Expc_all = {};
Expci_all = {};
% Start iterating until stopping criterion triggers
while true
    
    % Display iteration information
    if options.verbosity >= 2
        fprintf('%5d\t%9d\t%+.4e\t%.4e\n', iter, costevals, cost, gradnorm);
    end
    
    % Start timing this iteration
    timetic = tic();
    
    % Run standard stopping criterion checks
    [stop reason] = stoppingcriterion(problem, x, options, ...
        info, iter+1);
    % Run specific stopping criterion check
    if ~stop && stats.stepsize < options.minstepsize
        stop = true;
        reason = 'Last stepsize smaller than minimum allowed.';
    end
    
    if stop
        if options.verbosity >= 1
            fprintf([reason '\n']);
        end
        break;
    end
    
        
    if ~isfield(stats,'curiterratio')
        stats.curiterratio = 0.5;
    end

    if stats.curiterratio == 0 % This is the last iteration in a mini-batch
        [cost grad storedb] = getCostGrad(problem, x, storedb);
        gradnorm = problem.M.norm(x, grad);
        % Update BFGS inverse Hessian matrix and descent direction
        %  It is implemented by unrolling the inverse Hessian update
        if isempty(gd_all)
            desc_dir = problem.M.lincomb(x, 1/gradnorm, grad);
        elseif isfield(problem.M,'transpf')
            desc_dir = desc_dir_cal(grad, problem.M, grad_diff_all, ...
                desc_dir_all, x_all, ddgd_all, gd_all, length(gd_all) , H, ...
                Expc_all, Expci_all);
        else
            desc_dir = desc_dir_cal(grad, problem.M, grad_diff_all, ...
                desc_dir_all, x_all, ddgd_all, gd_all, length(gd_all) , H);
        end
        
        % Change search direction because it is gradient descend
        desc_dir = problem.M.lincomb(newx, -1, desc_dir);
    end
    
    % The line search algorithms require the directional derivative of the
    % cost at the current point x along the search direction.
    df0 = problem.M.inner(x, grad, desc_dir);
    if df0 > 0
        if options.verbosity >= 1
            fprintf(['Line search warning: got an ascent direction ' ...
                '(df0 = %2e), went the other way.\n'], df0);
        end
        desc_dir = problem.M.lincomb(x, -1, desc_dir);
        df0 = -df0;
    end
    
    % Execute line search
    [stepsize newx storedb lsmem lsstats] = options.linesearch(problem, ...
        x, desc_dir, cost, df0, options, storedb, lsmem, grad);

    % Compute the new cost-related quantities for x
    if isfield(lsmem,'grad')
        newcost = lsmem.cost;
        newgrad = lsmem.grad;
    else
        [newcost newgrad storedb] = getCostGrad(problem, newx, storedb);
    end
    
    costevals = costevals + lsstats.costevals;
    newgradnorm = problem.M.norm(newx, newgrad);
    
    % Make sure we don't use too much memory for the store database
    storedb = purgeStoredb(storedb, options.storedepth);
    
    % Using previous and new information to update
    % !!! This is important that parallel transport in this part should be
    %   associated with the retraction used in line-search (Ring&Wirth)
    if isfield(lsmem,'transg')
        gradC = lsmem.transg;
        %lsmem.newd = newd;
        [Expc, Expci] = problem.M.transpstore(x, newx);
    else
        gradC = problem.M.transp(x, newx, grad);
        if isfield(problem.M,'transpf')
            [Expc, Expci] = problem.M.transpstore(x, newx);
        end
    end
    %Expc.D{1}.sigmat
    grad_diff = problem.M.lincomb(newx, 1, newgrad, -1, gradC);
         
    
    % Parallel transport descent to the new point
    if isfield(lsmem,'newd')
        desc_dir_step = problem.M.lincomb(newx, stepsize, lsmem.newd);
        %problem.M.transpf(Expc, desc_dir_step);
    else
        % Multiplying the stepsize in descent direction
        desc_dir_step = problem.M.lincomb(x, stepsize, desc_dir);
        desc_dir_step = problem.M.transp(x, newx, desc_dir_step);
    end
    % disp('Calculating Updates');
    % Update the previous saved info
    if isfield(problem.M,'transpf')
        [grad_diff_all, desc_dir_all, x_all, gd_all, ddgd_all, H, ...
            Expc_all, Expci_all] = ...
            lbfgs_update(newx, problem.M, grad_diff, desc_dir_step, ...
            grad_diff_all, desc_dir_all, options.numgrad, x_all, gd_all, ...
            ddgd_all, H, options.verbosity, Expc, Expci, Expc_all, Expci_all);
    else
        [grad_diff_all, desc_dir_all, x_all, gd_all, ddgd_all, H] = ...
            lbfgs_update(newx, problem.M, grad_diff, desc_dir_step, ...
            grad_diff_all, desc_dir_all, options.numgrad, x_all, gd_all, ...
            ddgd_all, H, options.verbosity);
    end
    if H == 0
        break;
    end
    
    if stats.curiterratio < 1 % This is NOT the last iteration in a mini-batch
        % Update BFGS inverse Hessian matrix and descent direction
        %  It is implemented by unrolling the inverse Hessian update
        if isempty(gd_all)
            desc_dir = problem.M.lincomb(x, 1/newgradnorm, newgrad);
        elseif isfield(problem.M,'transpf')
            desc_dir = desc_dir_cal(newgrad, problem.M, grad_diff_all, ...
                desc_dir_all, x_all, ddgd_all, gd_all, length(gd_all) , H, ...
                Expc_all, Expci_all);
        else
            desc_dir = desc_dir_cal(newgrad, problem.M, grad_diff_all, ...
                desc_dir_all, x_all, ddgd_all, gd_all, length(gd_all) , H);
        end
        
        % Change search direction because it is gradient descend
        desc_dir = problem.M.lincomb(newx, -1, desc_dir);
        
    end
    % Update iterate info
    x = newx;
    cost = newcost;
    grad = newgrad;
    gradnorm = newgradnorm;
    
    % iter is the number of iterations we have accomplished.
    iter = iter + 1;
    % disp('Done');
    % Log statistics for freshly executed iteration
    stats = savestats();
    info(iter+1) = stats; %#ok<AGROW>
    
end

info = info(1:iter+1);

if options.verbosity >= 1
    fprintf('Total time is %f [s] (excludes statsfun)\n', ...
        info(end).time);
end



% Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = iter;
        stats.cost = cost;
        stats.gradnorm = gradnorm;
        if iter == 0
            stats.stepsize = nan;
            stats.time = toc(timetic);
            stats.linesearch = [];
        else
            stats.stepsize = stepsize;
            stats.time = info(iter).time + toc(timetic);
            stats.linesearch = lsstats;
        end
        stats = applyStatsfun(...
            problem, x, storedb, options, stats);
    end

end
