%% |mxe_options|
% *Note:* This is a private function.
%
% Generate default estimation options, Add given extra default options, and
% then set any given options.
%
% *Syntax*
%
%   options = mxe_options()
%   options = mxe_options(options)
%   options = mxe_options(options, extra_defaults)
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

function options = mxe_options(options, extra_defaults)

% Note: some options might not be applicable in some estimation functions.
% Note: changes should be reflected in documentation (estimation_options.m)

    if nargin < 2
        extra_defaults = [];
    end
        
    % don't fill the options if it is filled already and there are no extra options
    if nargin>0 && isfield(options, 'isfull__') && isempty(extra_defaults) 
        return
    end

    default_options = struct(...
        'isfull__', true, ... Flag used internally to prevent filling the default options more than once
        'solver', 'default', ... Solver handle, solver name or 'default' to call estimatedefault
        'theta0', [], ... Initial theta. If empty, the init function of the distribution is called to obtain an initial theta.
        'penalize', false, ... Add penalization to the cost function
        'penalizertheta', [], ... Penalizer parameters. If empty, uses the output of the penalizerparam function of the distribution on estimation data.
        'regularize', false, ... % Add regularization to the cost function.
        'minibatch', false, ... Flag to enable mini-batching, or a minibatch_options structure to customize mini-batching.
        'datapatchsize', Inf, ... Data is loaded patch-by-patch with patches of this size, when calculating the cost function. A value of Inf disables this feature.
        'checkgradient', false, ... Check gradient of the optimization problem instead of performing the full estimation
        'verbosity', 1, ... Controls how much information is displayed during estimation (0: none; 1: at init and exit; 2: at every iteration)
        'previnfo', [], ... The info structure array of a previous estimation function to continue from.
        ... Plotting
        'plotcost', false, ... Flag to enable plotting the cost value during estimation, or a plot_options structure to customize the plot
        'plotgradnorm', false, ... Flag to enable plotting the norm of the gradient during estimation, or a plot_options structure to customize the plot
        'visualization', false, ... Flag to enable visualization of the estimation process for 2-D or 3-D data
        ... Stopping Criteria
        'crossval', false, ... Flag to enable cross-validation, or a crossval_options structure to customize cross-validation.
        'miniter', 0, ... Minimum number of iterations
        'maxiter', Inf, ... Maximum number of iterations
        'maxtime', Inf, ... Maximum running time, in seconds
        'maxcostevals', Inf, ... Maximum number of evaluations of the cost function
        'tolcost', -Inf, ... Stop as soon as the cost drops below this tolerance
        'tolgradnorm', 1e-6, ... Stop as soon as the norm of the gradient drops below this tolerance
        'minstepsize', 1e-10, ... Stop as soon as the norm of the step size drops below this tolerance
        'tolcostdiff', 1e-6, ... Stop as soon as the decrease in cost in an iteration is less than this tolerance
        ... Debug
        'debug', false, ... Set to true to display additional information for debugging purposes.
        ... Split and Merge
        'sm', struct(...
            'splitcriterion', 'kl', ... 'kl', 'mll', 'entropy', 'rand'
            'mergecriterion', 'kl', ... 'kl', 'overlap', 'rand'
            'splitinit', 'default', ... Split method name or function handle or e.g. options.splitinit.mvn.mu = 'method'
            'mergeinit', 'default' ... Split method name or function handle or e.g. options.mergeinit.mvn.mu = 'method'
            ), ...
        ... stochastic gradient descend
        'sgd', struct(...
            'batchnum', 100, ... Number of batches in SGD
            'epoch', 50, ... Stop as soon as number of epochs reachs epoc
            'stepsize', 0.1 ... step size in SGD
            ) ...
        ... See also mxe_inneroptions
        );
    
    
    % add extra default options
    if ~isempty(extra_defaults)
        default_options = mxe_setfields(default_options, extra_defaults);
    end
    
    
    % set the given options
    if nargin < 1
        options = default_options;
    else
        % theta0 should not become lower case
        if isfield(options, 'theta0')
            theta0 = options.theta0;
        else
            theta0 = [];
        end
        options = mxe_setfields(default_options, options);
        options.theta0 = theta0;
    end
    
    % convert the options to standard form
    options = normalize_options(options);
end



function options = normalize_options(options)

    options = normalize_string_options(options);
    options = normalize_minibatch_options(options);
    options = normalize_crossval_options(options);
    options = normalize_plot_options(options);
    options = normalize_visualization_options(options);
end



function options = normalize_string_options(options)
% convert given string options to standard form

    if isfield(options, 'solver') && ischar(options.solver)
        options.solver = lower(options.solver);
    end
    if isfield(options, 'splitcriterion') && ischar(options.splitcriterion)
        options.splitcriterion = lower(options.splitcriterion);
    end
    if isfield(options, 'mergecriterion') && ischar(options.mergecriterion)
        options.mergecriterion = lower(options.mergecriterion);
    end
    if isfield(options, 'splitinit') && ischar(options.splitinit)
        options.splitinit = lower(options.splitinit);
    end
    if isfield(options, 'mergeinit') && ischar(options.mergeinit)
        options.mergeinit = lower(options.mergeinit);
    end
end


    
function options = normalize_minibatch_options(options)
% make the minibatch options standard (full struct)

    default_minibatch_options = struct( ...
        'size', 1000, ... Number of data points in each minibatch. A value of zero disables mini-batching
        'iter', 10, ... Number of iterations per mini-batch
        'discardhistory', true, ... Discard optimization history and start again from the final theta before moving to the next mini-batch
        'overlap', 0 ... Number of data points shared between two consecutive mini-batches
        );

    if islogical(options.minibatch)
        enabled = options.minibatch;
        options.minibatch = default_minibatch_options;
        if ~enabled
            options.minibatch.size = 0;
        end
    else
        options.minibatch = mxe_setfields(default_minibatch_options, options.minibatch);
    end
end


    
function options = normalize_crossval_options(options)
% make the cross-validation options standard (full struct)

    default_crossval_options = struct( ...
        'enabled', true, ... Flag to enable cross-validation
        'fraction', 0.2, ... Fraction of data to be used for cross validation. The rest of data is used for training.
        'toliter', 10 ... Permitted number of iterations to have a rising cross-validation cost
        );

    if islogical(options.crossval)
        enabled = options.crossval;
        options.crossval = default_crossval_options;
        options.crossval.enabled = enabled;
    else
        options.crossval = mxe_setfields(default_crossval_options, options.crossval);
    end
end


    
function options = normalize_plot_options(options)
% make the plotting options standard (full struct)
    
    default_plot_options = struct( ...
        'enabled', true, ... Flag to enable the plot
        'axes', 0, ... Handle of the axes object for the plot. A value of zero means creating a new figure and axes.
        'avgiter', 1, ... Use averaged values over this number of iterations for plotting
        'itercount', 100, ... Plot over this number of last iterations
        'xlabel', 'iterations', ... X-axis label
        'ylabel', '', ... Y-axis label
        'title', '', ... Plot title
        'legend', 'default', ... Cell array of labels for the legend, or 'default' to use the default labels
        'log', true, ... Flag to use a logarithmic Y-scale
        'grid', true ... Flag to show grid lines
        );
    
    % cost plot
    if islogical(options.plotcost)
        enabled = options.plotcost;
        options.plotcost = default_plot_options;
        options.plotcost.enabled = enabled;
    else
        options.plotcost = mxe_setfields(default_plot_options, options.plotcost);
    end
    if options.plotcost.enabled
        options.plotcost.axes = check_valid_handle(options.plotcost.axes, 'Cost');
    end

    % gradient norm plot
    if islogical(options.plotgradnorm)
        enabled = options.plotgradnorm;
        options.plotgradnorm = default_plot_options;
        options.plotgradnorm.enabled = enabled;
    else
        options.plotgradnorm = mxe_setfields(default_plot_options, options.plotgradnorm);
    end
    if options.plotgradnorm.enabled
        options.plotgradnorm.axes = check_valid_handle(options.plotgradnorm.axes, 'Gradient Norm');
    end

end



function options = normalize_visualization_options(options)
% make the visualization options standard (full struct)
    
    default_visualization_options = struct( ...
        'enabled', true, ... Flag to enable visualization
        'axes', 0, ... Handle of the axes object for the visualization. A value of zero means creating a new figure and axes.
        'stoponclose', false, ... Flag indicating if we should stop estimation on closing the visualization figure
        'dataplottype', 'patch', ... Use 'patch' plot (faster) or 'scatter' plot (slower) for data visualization
        'visobject', [] ... A visualization object (returned from mxe_visualization) to use instead of creating a new one
        );
    
    if islogical(options.visualization)
        enabled = options.visualization;
        options.visualization = default_visualization_options;
        options.visualization.enabled = enabled;
    else
        options.visualization = mxe_setfields(default_visualization_options, options.visualization);
    end
    if options.visualization.enabled
        options.visualization.axes = check_valid_handle(options.visualization.axes, 'Visualization');
    end
end



function ax = check_valid_handle(ax, figure_name)

    if ~( ishandle(ax) && strcmp(get(ax,'type'),'axes') )
        figure('Name', figure_name)
        ax = gca;
    end
end
