function [stepsize newx storedb lsmem lsstats] = ...
       linesearch_wolfe(problem, x, d, f0, df0, options, storedb, lsmem)
% Adaptive line search algorithm (step size selection) for descent methods.
%
% function [stepsize newx storedb lsmem lsstats] = 
%      linesearch_adaptive(problem, x, d, f0, df0, options, storedb, lsmem)
%
% linesearch algorithm satisfying strong wolfe conditions.
% See Equations 1a and 2 of 
% W. Ring, B. Wirth, " Optimization methods on Riemannian manifolds and
%         their application to shape space"
%    x_new = retr(x,d,alpha)
%    d_new = transp(x,x_new,d)
%    f_new < f0 + alpha*suff_decr*df0 (Armijo Condiiton) 
%    abs(inner(x_new,df_new,d_new)) <= -wolfe_cte* inner(x,df0,d) 
%                                               (Strong Wolfe Condition)
% 
% Based on the algorithm explained in chapter 2 of 
%   "Numerical Optimization" by Joerg Nocedal & Martin Wainwright
%
% Inputs/Outputs : see help for linesearch
%
% See also: steepestdescent conjugategradients linesearch

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Original author: Reshad Hosseini, May. 19, 2013.
%
% Inpired by the linesearch code of minFunc toolbox


   
        
    % Wolfe default parameters. These can be overwritten in the
    % options structure which is passed to the solver.
    default_options.ls_optimism = 1.1;
    default_options.ls_init_step_proc = 1;
    default_options.ls_suff_decr = 1e-4;
    default_options.ls_wolfe_cte = 0.9;
    default_options.ls_max_steps = 20;
    default_options.ls_initial_stepsize = 1;
    
    % When extrapolation leads to inf log-likelihood applies Armijo rule
    % with certain contraction factor
    default_options.ls_contraction_factor = .1;
    options = mergeOptions(default_options, options);
    
    suff_decr = options.ls_suff_decr;
    init_step_proc = options.ls_init_step_proc; 
    max_ls_steps = options.ls_max_steps;
    wolfe_cte = options.ls_wolfe_cte;
    contraction_factor = default_options.ls_contraction_factor;
    optimism = options.ls_optimism;
    new_alpha = options.ls_initial_stepsize;
    M = problem.M;
    
    if init_step_proc == 1
        if isfield(lsmem, 'f0')
            % Pick initial step size based on where we were last time,
            new_alpha = 2*(f0-lsmem.f0)/ df0;
            % In order to Gaurantee superlinear convergence.
            new_alpha = min(1, optimism * new_alpha);
        end
    end
    
    if init_step_proc == 2
        if isfield(lsmem, 'df0')
            % Pick initial step size based on where we were last time,
            new_alpha = lsmem.alpha * lsmem.df0 / df0;
            new_alpha = min(1, optimism * new_alpha);
        end
    end
    % Bracketing Phase
    ls_iter = 0;
    done = false;
    do_plot = false;
    alpha = 0;
    debug = false;
    f = f0;
    newx = problem.M.retr(x, d, new_alpha);
    first_run = false; %Do not Evaluate Gradient in first runs
    if first_run
        new_df = -1;
        [new_f storedb] = getCost(problem, newx, storedb);
    else
        [new_f new_grad storedb] = getCostGrad(problem, newx, storedb);
        if ~is_illegal(new_f)
            new_df = calcDf(M, x, newx, d, new_grad);
        end
    end

    if ~is_illegal(new_f)
        points = [alpha f df0; new_alpha new_f new_df];
    end
    
    if debug
        disp('start...');
    end
    bestf = -Inf;
    while ls_iter < max_ls_steps
        ls_iter = ls_iter + 1;
        % Reaching a point with Illegal function value -> Do Armijo Rule
        if is_illegal(new_f)
            if debug
                disp('illegal');
            end
            while true
                % Reduce the step size,
                new_alpha = contraction_factor * new_alpha;
                % and look closer down the line
                newx = problem.M.retr(x, d, new_alpha);
                [new_f storedb] = getCost(problem, newx, storedb);
                ls_iter = ls_iter + 1;
                if new_f <= f0 + suff_decr*new_alpha*df0
                    [new_grad storedb] = getGradient(problem, newx, storedb);
                    if ~is_illegal(new_grad)
                        break;
                    end
                end
                % Make sure we don't run out of budget
                if ls_iter >= max_ls_steps 
                    disp('beyond number of steps');
                    if is_illegal(new_f) || new_f > f0
                        newx = x;
                        [new_f new_grad storedb] = getCostGrad(problem, newx, storedb);
                    else
                        [new_grad storedb] = getGradient(problem, newx, storedb);
                    end
                    break;
                end
            end
            done = true;
            break;
        end
        % Different conditions and choosing bracket accordingly
        if new_f > f0 + new_alpha*suff_decr*df0 || (ls_iter > 1 && new_f >= f)
            if debug
                disp('1st');
            end
            break;
        else
            if debug
                disp('No 1st');
            end
            if first_run
                [new_grad storedb] = getGradient(problem, newx, storedb);
                new_df = calcDf(M, x, newx, d, new_grad);
                points = [alpha f df0; new_alpha new_f new_df];
                first_run = false;
            end
            if abs(new_df) <= -wolfe_cte*df0
                if debug
                    disp('2nd');
                end
                done = true;
                break;
            elseif new_df >= 0
                if debug
                    disp('3rd');
                end
                break;
            end
        end 
        if new_f < f0 + new_alpha*suff_decr*df0
            bestf = new_f;
        end
        % Cubic Extrapolation        
        min_alpha = new_alpha + 0.01*(new_alpha-alpha);
        max_alpha = new_alpha * 10;
        %points = [alpha f df; new_alpha new_g new_df];
        alpha = new_alpha;
        new_alpha = polyinterp(points, do_plot, min_alpha, max_alpha);
        % moving along geodesic to find the new point
        newx = problem.M.retr(x, d, new_alpha);
        f = new_f;
        df = new_df;
        [new_f new_grad storedb] = getCostGrad(problem, newx, storedb);
        new_df = calcDf(M, x, newx, d, new_grad);
        points = [alpha f df; new_alpha new_f new_df];
    end
    
    % use index for the place of higher value and lower value
    if new_f > f
        index_low = 1;
        index_high = 2;
        low_f = f;
    else
        index_low = 2;
        index_high = 1;
        low_f = new_f;
    end
    
    % Zooming Phase
    while ~done && ls_iter < max_ls_steps
        ls_iter = ls_iter + 1;
        % Cubic Extrapolation
        new_alpha = polyinterp(points, do_plot);
        
        % Test that we are making sufficient progress
        interval = max(points(:,1))-min(points(:,1));
        if min(max(points(:,1))-new_alpha,new_alpha-min(points(:,1)))...
                /interval < 0.1
            if debug
                fprintf('Interpolation close to boundary');
                fprintf(', Evaluating at 0.1 away from boundary\n');
            end
            if abs(new_alpha-max(points(:,1))) < abs(new_alpha-min(points(:,1)))
                new_alpha = max(points(:,1)) - 0.1*interval;
            else
                new_alpha = min(points(:,1)) + 0.1*interval;
            end
        end
        newx = problem.M.retr(x, d, new_alpha);
        if first_run
            [new_f storedb] = getCost(problem, newx, storedb);
            new_df = -1;
        else
            [new_f new_grad storedb] = getCostGrad(problem, newx, storedb);
            new_d = M.transp(x, newx, d);
            new_df = M.inner(newx, new_d, new_grad);
        end

        if (new_f > f0 + new_alpha*suff_decr*df0 || new_f >= low_f) && new_f > bestf
            if debug
                disp('4th');
            end
            points(index_high, :) = [new_alpha new_f new_df];
            % put new points in the situation of high point
        else
            if first_run
                [new_grad storedb] = getGradient(problem, newx, storedb);
                new_d = M.transp(x, newx, d);
                new_df = M.inner(newx, new_d, new_grad);
                first_run = false;
            end
            if abs(new_df) <= -wolfe_cte*df0
                if debug
                    disp('5th');
                end
                % Wolfe conditions satisfied
                done = true;
            elseif new_df*(points(index_high,1)-points(index_low,1)) >= 0
                % new point is the low point
                if debug
                    disp('6th');
                end
                points(index_high, :) = points(index_low, :);
            else
                if debug
                    disp('7th');
                end
                %points(index_low, :) = [new_alpha new_f new_df];
            end
            points(index_low, :) = [new_alpha new_f new_df];
            low_f = new_f;
        end
            
        if abs((points(1,1)-points(2,1))*new_df) < 1e-10 && new_f < f0
            if debug
                fprintf('Line Search can not make further progress\n');
            end
            break;
        end
    end
    
    % Return the cost and gradient in the new point
    lsmem.cost = new_f;
    lsmem.grad = new_grad;
    lsmem.f0 = f0;
    lsmem.df0 = df0;
    lsmem.alpha = new_alpha;
    if debug
        ls_iter
    end
    % As seen outside this function, stepsize is the size of the vector we
    % retract to make the step from x to newx. 
    stepsize = new_alpha;
    
    % Save some statistics also, for possible analysis.
    lsstats.costevals = ls_iter;
    lsstats.stepsize = stepsize;
    lsstats.alpha = new_alpha;
    
end

function ilegal = is_illegal(x)
    x = obj2vec(x);
    ilegal = any(imag(x(:))) | any(isnan(x(:))) | any(isinf(x(:)));
end

function  new_df = calcDf(M, x, newx, d, new_grad)
    new_d = M.transp(x, newx, d);
    new_df = M.inner(newx, new_d, new_grad);
end
