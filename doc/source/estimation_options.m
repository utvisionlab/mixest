%% Estimation Options
% All estimation functions in MixEst accept an optional |options| structure
% to let the user control the behavior of the estimation algorithm. Each
% field of the structure corresponds to an option. The default value is
% used for any missing option.
%
% The estimation functions read only the fields that they understand from
% the structure and leave the other fields intact, possibly passing them to
% other functions that are called inside (such as Manopt solvers). The list
% of options available for each estimation function is mentioned in its
% documentation.
%

%% List of commonly-used options
% Here is a list of common options used in estimation functions. For
% information about the available options for a specific function, refer to
% the documentation of that function. Note that the option names are NOT
% case sensitive and are all converted to lower case before processing.
%

clear
options = mxe_options();

k = 1;
T(k).Name = 'solver';
T(k).Type = 'string / function handle';
T(k).Default_Value = ['''' options.(lower(T(k).Name)) ''''];
T(k).Description = 'Solver handle, solver name or ''default'' to call estimatedefault. See the <a href="#2">next section</a> for more details.';

k = k + 1;
T(k).Name = 'theta0';
T(k).Type = 'parameter structure';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'Initial theta. If empty, the init function of the distribution is called to obtain an initial theta';

k = k + 1;
T(k).Name = 'penalize';
T(k).Type = 'logical';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'Add penalization to the cost function';

k = k + 1;
T(k).Name = 'penalizerTheta';
T(k).Type = 'parameter structure';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'Penalizer parameters. If empty, uses the output of the penalizerparam function of the distribution on estimation data.';

k = k + 1;
T(k).Name = 'miniBatch';
T(k).Type = 'logical / structure';
T(k).Default_Value = false; %options.(lower(T(k).Name));
T(k).Description = 'Flag to enable mini-batching, or a minibatch_options structure to customize mini-batching. See <a href="#3">below</a> for more details.';

k = k + 1;
T(k).Name = 'dataPatchSize';
T(k).Type = 'integer';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'Data is loaded patch-by-patch with patches of this size, when calculating the cost function. A value of Inf disables this feature.';

k = k + 1;
T(k).Name = 'checkGradient';
T(k).Type = 'logical';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'Check gradient of the optimization problem instead of performing the full estimation';

k = k + 1;
T(k).Name = 'verbosity';
T(k).Type = 'integer';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'Controls how much information is displayed during estimation (0: none; 1: at init and exit; 2: at every iteration)';

k = k + 1;
T(k).Name = 'prevInfo';
T(k).Type = 'structure array';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'The info structure array of a previous estimation function to continue from.';

k = k + 1;
T(k).Name = 'plotCost';
T(k).Type = 'logical / structure';
T(k).Default_Value = false; %options.(lower(T(k).Name));
T(k).Description = 'Flag to enable plotting the cost value during estimation, or a plot_options structure to customize the plot. See <a href="#4">below</a> for more details.';

k = k + 1;
T(k).Name = 'plotGradNorm';
T(k).Type = 'logical / structure';
T(k).Default_Value = false; %options.(lower(T(k).Name));
T(k).Description = 'Flag to enable plotting the norm of the gradient during estimation, or a plot_options structure to customize the plot. See <a href="#4">below</a> for more details.';

k = k + 1;
T(k).Name = 'visualization';
T(k).Type = 'logical / structure';
T(k).Default_Value = false; %options.(lower(T(k).Name));
T(k).Description = 'Flag to enable 2-D or 3-D visualization during estimation, or a vis_options structure to customize visualization. See <a href="#5">below</a> for more details.';

k = k + 1;
T(k).Name = 'crossVal';
T(k).Type = 'logical / structure';
T(k).Default_Value = false; %options.(lower(T(k).Name));
T(k).Description = 'Flag to enable cross-validation, or a crossval_options structure to customize cross-validation. See <a href="#6">below</a> for more details.';

k = k + 1;
T(k).Name = 'minIter';
T(k).Type = 'integer';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'Minimum number of iterations';

k = k + 1;
T(k).Name = 'maxIter';
T(k).Type = 'integer';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'Maximum number of iterations';

k = k + 1;
T(k).Name = 'maxTime';
T(k).Type = 'double';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'Maximum running time, in seconds';

k = k + 1;
T(k).Name = 'maxCostEvals';
T(k).Type = 'integer';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'Maximum number of evaluations of the cost function';

k = k + 1;
T(k).Name = 'tolCost';
T(k).Type = 'double';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'Stop as soon as the cost drops below this tolerance';

k = k + 1;
T(k).Name = 'tolGradNorm';
T(k).Type = 'double';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'stop as soon as the norm of the gradient drops below this tolerance';

k = k + 1;
T(k).Name = 'minStepSize';
T(k).Type = 'double';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'Stop as soon as the norm of the step size drops below this tolerance';

k = k + 1;
T(k).Name = 'tolCostDiff';
T(k).Type = 'double';
T(k).Default_Value = options.(lower(T(k).Name));
T(k).Description = 'Stop as soon as the decrease in cost in an iteration is less than this tolerance';

k = k + 1;
T(k).Name = 'statsfun';
T(k).Type = 'function handle';
T(k).Default_Value = '-';
T(k).Description = 'Used to log custom information during estimation. See <a href="#7">below</a> for more details.';

k = k + 1;
T(k).Name = 'stopfun';
T(k).Type = 'function handle';
T(k).Default_Value = '-';
T(k).Description = 'Used to define custom termination criteria for the estimation procedure. See <a href="#8">below</a> for more details.';

k = k + 1;
T(k).Name = 'costgrad';
T(k).Type = 'function handle';
T(k).Default_Value = '-';
T(k).Description = 'Used to define custom cost function for the estimation procedure. See <a href="#9">below</a> for more details.';

disp(mxe_htmltable(T))

%% Solvers
% In the main estimation function, we use <http://www.manopt.org Manopt>
% solvers to solve the optimization problem involved in model parameter
% estimation. You can select the employed solver by setting
% |options.solver| to the handle of the solver function (e.g.
% |options.solver = @trustregions;|), or alternatively you can set
% |options.solver| to one of the string values from the following table.
%

clear T

k = 1;
T(k).Value = '''default''';
T(k).Description = 'If the distribution structure exposes an estimatedefault function, calls it. Otherwise uses the LBFGS solver';

k = k + 1;
T(k).Value = '''lbfgs''';
T(k).Description = 'LBFGS solver';

k = k + 1;
T(k).Value = '''cg'' or ''conjugategradient''';
T(k).Description = 'conjugate-gradient solver';

k = k + 1;
T(k).Value = '''tr'' or ''trustregions''';
T(k).Description = 'trust-regions solver';

k = k + 1;
T(k).Value = '''sd'', ''gd'', ''steepestdescent'' or ''gradientdescent''';
T(k).Description = 'steepest-descent (gradient-descent) solver';

disp(mxe_htmltable(T))

%% Mini-batching options
% To enable mini-batching, you can set |options.miniBatch| to |true| to use
% the default mini-batching options. Alternatively you may set special
% mini-batching options as fields in |options.miniBatch| (as a structure)
% to customize mini-batching. For example the following code sets the
% number of data points in each mini-batch to 400 and the number of
% iterations per mini-batch to 20. The other mini-batcing options will be
% set to their default values.
%
%   options.minibatch.size = 400;
%   options.minibatch.iter = 20;
%
% The following table lists the mini-batching options:
%

clear T

k = 1;
T(k).Name = 'size';
T(k).Type = 'integer';
T(k).Default_Value = 1000; %options.minibatch.(lower(T(k).Name));
T(k).Description = 'Number of data points in each minibatch. A value of zero disables mini-batching';

k = k + 1;
T(k).Name = 'iter';
T(k).Type = 'integer';
T(k).Default_Value = options.minibatch.(lower(T(k).Name));
T(k).Description = 'Number of iterations per mini-batch';

k = k + 1;
T(k).Name = 'discardHistory';
T(k).Type = 'logical';
T(k).Default_Value = options.minibatch.(lower(T(k).Name));
T(k).Description = 'Discard optimization history and start again from the final theta before moving to the next mini-batch';

k = k + 1;
T(k).Name = 'overlap';
T(k).Type = 'integer';
T(k).Default_Value = options.minibatch.(lower(T(k).Name));
T(k).Description = 'Number of data points shared between two consecutive mini-batches';

disp(mxe_htmltable(T))

%% Plotting options
% MixEst gives you the option to plot the optimization cost and/or gradient
% norm during estimations. You can set |options.plotCost| or
% |options.plotGradNorm| to |true| to enable this feature. Alternatively
% you may set predefined field(s) in these options (as structures) to
% customize plotting behavior. The following example shows how to set
% options to show the cost and gradient norm plot in a single figure during
% estimation.
%
%   figure
%   options.plotcost.axes = subplot(2,1,1);
%   options.plotgradnorm.axes = subplot(2,1,2);
%
% The following table lists the plotting options:
%

clear T

k = 1;
T(k).Name = 'enabled';
T(k).Type = 'logical';
T(k).Default_Value = true; %options.plotcost.(lower(T(k).Name));
T(k).Description = 'Flag to enable the plot';

k = k + 1;
T(k).Name = 'axes';
T(k).Type = 'axes handle';
T(k).Default_Value = options.plotcost.(lower(T(k).Name));
T(k).Description = 'Handle of the axes object for the plot. A value of zero means creating a new figure and axes.';

k = k + 1;
T(k).Name = 'avgiter';
T(k).Type = 'integer';
T(k).Default_Value = options.plotcost.(lower(T(k).Name));
T(k).Description = 'Use averaged values over this number of iterations for plotting';

k = k + 1;
T(k).Name = 'iterCount';
T(k).Type = 'integer';
T(k).Default_Value = options.plotcost.(lower(T(k).Name));
T(k).Description = 'Plot over this number of last iterations';

k = k + 1;
T(k).Name = 'xlabel';
T(k).Type = 'string';
T(k).Default_Value = ['''' options.plotcost.(lower(T(k).Name)) ''''];
T(k).Description = 'X-axis label';

k = k + 1;
T(k).Name = 'ylabel';
T(k).Type = 'string';
T(k).Default_Value = ['''' options.plotcost.(lower(T(k).Name)) ''''];
T(k).Description = 'Y-axis label';

k = k + 1;
T(k).Name = 'title';
T(k).Type = 'string';
T(k).Default_Value = ['''' options.plotcost.(lower(T(k).Name)) ''''];
T(k).Description = 'Plot title';

k = k + 1;
T(k).Name = 'legend';
T(k).Type = 'cell array of strings';
T(k).Default_Value = ['''' options.plotcost.(lower(T(k).Name)) ''''];
T(k).Description = 'Cell array of labels for the legend, or ''default'' to use the default labels';

k = k + 1;
T(k).Name = 'log';
T(k).Type = 'logical';
T(k).Default_Value = options.plotcost.(lower(T(k).Name));
T(k).Description = 'Flag to use a logarithmic Y-scale';

k = k + 1;
T(k).Name = 'grid';
T(k).Type = 'logical';
T(k).Default_Value = options.plotcost.(lower(T(k).Name));
T(k).Description = 'Flag to show grid lines';

disp(mxe_htmltable(T))

%% Visualization options
% We have implemented simple visualizations for some of the distributions
% that might be useful for you (They were really helpful for us during
% development). If you have 2-D or 3-D data, You can set
% |options.visualization| to |true| to give it a try. Alternatively you may
% set predefined field(s) in |options.visualization| (as a structure) to
% customize the visualization. This feature is experimental and is not
% fully documented yet but the following example illustrates how you can
% set some visualization options.
%
%   % main visualization options
%   options.visualization.axes = gca;
%   options.visualization.stopOnClose = true; % stop estimation when the figure is closed
%   options.visualization.dataPlotType = 'patch'; % or 'scatter'
%   % visualization options for the mixture distribution
%   options.visualization.mixture.colorize = true;
%   options.visualization.mixture.colorMap = hsv(3);
%   options.visualization.mixture.showLabels = true;
%   % visualization options for the MVN distribution
%   options.visualization.mvn.showCenter = false;
%   % load some data and perform estimation
%   load data2d
%   D = mixturefactory(mvnfactory(2), 3);
%   D.estimate(data, options)
%

%% Cross-validation options
% To enable cross-validation, you can set |options.crossVal| to |true| to
% use the default cross-validation options. Alternatively you may set
% special cross-validation options as fields in |options.crossVal| (as a
% structure) to customize cross-validation. For example the following code
% asks for half of the data to be used for cross-validation. The other
% cross-validation options will be set to their default values.
%
%   options.crossval.fraction = 0.5;
%
% Note that your data should be shuffled before using cross-validation.
%
% The following table lists the cross-validation options:
%

clear T

k = 1;
T(k).Name = 'enabled';
T(k).Type = 'logical';
T(k).Default_Value = true; %options.crossval.(lower(T(k).Name));
T(k).Description = 'Flag to enable cross-validation';

k = k + 1;
T(k).Name = 'fraction';
T(k).Type = 'double';
T(k).Default_Value = options.crossval.(lower(T(k).Name));
T(k).Description = 'Fraction of data to be used for cross validation. The rest of data is used for training.';

k = k + 1;
T(k).Name = 'tolIter';
T(k).Type = 'integer';
T(k).Default_Value = options.crossval.(lower(T(k).Name));
T(k).Description = 'Permitted number of iterations to have a rising cross-validation cost';

disp(mxe_htmltable(T))

%% |statsfun|
% If you add a |statsfun| field containing a function handle with prototype
% |stats = statsfun(D, theta, stats)| to the |options| structure, it will
% be called after each iteration completes with the current distribution
% structure |D| and its parameter values |theta| and the statistics
% structure |stats| that will be logged in the |info| structure array as
% described <estimation_statistics_structure.html here>, and gives you a
% chance to modify the stats structure. (This is based on how
% <http://www.manopt.org Manopt> declares |statsfun| but omitting its
% |problem| structure).
%
% *Example*
%
% Use the following code to log all the points (parameter values) visited
% during the estimation process:
%
%   options.statsfun = @mystatsfun;
%   function stats = mystatsfun(theta, stats)
%       stats.theta = theta;
%   end
%
% This will log all the points (parameter values) visited during the
% estimation process in the |info| struct-array returned by the estimation
% function.
%

%% |stopfun|
% If you add a |stopfun| field containing a function handle with prototype
% |stopnow = stopfun(D, theta, info, last)| to the |options| structure, it
% will be called after each iteration completes with the current
% distribution structure |D| and its parameter values |theta|, the whole
% |info| structure array built so far and an index |last| such that
% |info(last)| is the statistics structure corresponding to the current
% iteration (this is because |info| is pre-allocated, so that info(end)
% typically does not refer to the current iteration). The returned value
% |stopnow| should be a logical. If it is returned as |true|, the algorithm
% will terminate. (This is based on how <http://www.manopt.org Manopt>
% declares |stopfun| but omitting its |problem| structure).
%
% *Example*
%
%   options.stopfun = @mystopfun;
%   function stopnow = mystopfun(D, theta, info, last)
%       stopnow = (last >= 3 && info(last-2).cost - info(last).cost < 1e-3);
%   end
%
% This will tell the algorithm to exit as soon as two successive iterations
% combined have decreased the cost (negative log-likelihood) by less than
% |1e-3|.
%

%% |costgrad|
% If you add a |costgrad| field containing a function handle with prototype
% |[cost, grad] = costgrad(theta)| to the |options| structure, it will be
% used instead of the default cost/grad function (|mxe_costgrad.m|) to
% calculate the cost |cost| and/or its Riemannian gradient |grad| at every
% point |theta| during the estimation.
%
% *Example*
%
% The following example overrides the default cost/grad function and
% implements a preliminary one:
%
%   function test
%       D = mvnfactory(1);
%       data = randn(1, 1000);
%       options.solver = 'cg';
%       options.theta0 = D.randparam();
%       options.costgrad = @mycostgrad;
%       D.estimate(data, options);
%       %
%       %
%       function [cost, grad] = mycostgrad(theta)
%           [ll, store] = D.ll(theta, data);
%           cost = -ll;
%           if nargout > 1
%               dll = D.llgrad(theta, data, store); 
%               grad = D.M.egrad2rgrad(theta, dll); % convert the gradient to Riemannian
%               grad = D.M.lincomb(theta, -1, grad); % negate the gradient
%           end
%       end
%   end
%

%% Split-and-merge options
% Options for split-and-merge algorithms can be set through fields in the
% |options.sm| structure. The following options are common. Additional
% options for each algorithm are described in the documentation of the
% respective function.
%
% * *|splitCriterion|* (default |'kl'|) : Criterion for selecting a
% component to be split. Following values are available:
%
% <html>
%     <div class="table-responsive">
%         <table class="table table-striped table-bordered">
%             <thead>
%                 <tr>
%                     <th>Value</th>
%                     <th>Description</th>
%                 </tr>
%             </thead>
%             <tbody>
%                 <tr>
%                     <td><tt>'kl'</tt></td>
%                     <td>Maximum local KL-divergence between component distribution and empirical distribution</td>
%                 </tr>
%                 <tr>
%                     <td><tt>'mll'</tt></td>
%                     <td>Minimum mean local log-likelihood</td>
%                 </tr>
%                 <tr>
%                     <td><tt>'entropy'</tt></td>
%                     <td>Maximum entropy</td>
%                 </tr>
%                 <tr>
%                     <td><tt>'rand'</tt></td>
%                     <td>Random split</td>
%                 </tr>
%             </tbody>
%         </table>
%     </div>
% </html>
%
% * *|mergeCriterion|* (default |'kl'|) : Criterion for selecting two
% components to be merged. Following values are available:
%
% <html>
%     <div class="table-responsive">
%         <table class="table table-striped table-bordered">
%             <thead>
%                 <tr>
%                     <th>Value</th>
%                     <th>Description</th>
%                 </tr>
%             </thead>
%             <tbody>
%                 <tr>
%                     <td><tt>'kl'</tt></td>
%                     <td>Minimum symmetric KL-divergence</td>
%                 </tr>
%                 <tr>
%                     <td><tt>'overlap'</tt></td>
%                     <td>Maximum distribution overlap</td>
%                 </tr>
%                 <tr>
%                     <td><tt>'rand'</tt></td>
%                     <td>Random merge</td>
%                 </tr>
%             </tbody>
%         </table>
%     </div>
% </html>
%
% * *|splitInit|* (default |'default'|) : Initialization method for the
% components resulting from a split. You can set this field to 'default' to
% use the default methods for each parameter, or you may give a function
% handle implementing your custom initialization. The referred function
% should have the following prototype:
%
%   function [newtheta, store] = splitinit(D, idx, theta, options, data, store)
%
% where |D| is the mixture distribution structure, |idx| is the index of
% the component to be split, |theta| is the parameter structure for the
% mixture distribution, |options| is the full options structure, |data| is
% the <data_input.html given data>, and |store| can be used for
% <caching.html caching purposes>. The output |newtheta| is the parameter
% structure of the new mixture after the split (containing the parameters
% for an additional component). You should put the initialized parameter
% values for the first splitted component at |newtheta.D{idx}| and for the
% second in |newtheta.D{end}|. The component weights in |newtheta.p| should
% also be updated accordingly.
%
% * *|mergeInit|* (default |'default'|) : Initialization method for the
% component resulting from a merge. You can set this field to 'default' to
% use the default methods for each parameter, or you may give a function
% handle implementing your custom initialization. The referred function
% should have the following prototype:
%
%   function [newtheta, store] = mergeinit(D, idx1, idx2, theta, options, data, store)
%
% where |D| is the mixture distribution structure, |idx1| and |idx2| are
% the indices of the components to be merged, |theta| is the parameter
% structure for the mixture distribution, |options| is the full options
% structure, |data| is the <data_input.html given data>, and |store| can be
% used for <caching.html caching purposes>. The output |newtheta| is the
% parameter structure of the new mixture after the merge (with one less
% item than |theta|). You should put the initialized parameter values for
% the merged component at |newtheta.D{idx1}| and remove |theta.D{idx2}| in
% |newtheta|. The component weights in |newtheta.p| should also be updated
% accordingly.
%
% *Example*
%
% If you want the split-and-merge algorithm to select the component with
% maximum entropy when it requires to split a component, use the following
% option:
%
%   options.sm.splitCriterion = 'entropy';
%
% *Example*
%
% The following code shows how to use |options.sm.splitInit| to manually
% initialize the parameters of splitted components for a mixture of
% Gaussians. For the mean vectors, it uses k-means clustering (with k=2) on
% 10 data points that have the highest posterior probability under the
% original component. Each covariance matrix is initialized as a unit
% matrix with the same volume as the original component. The weight for
% each of the splitted components is set to half the weight of the original
% component.
%
%   options.sm.splitinit = @splitinit;
%   function [newtheta, store] = splitinit(D, idx, theta, options, data, store)
%     %
%     % add parameters for an additional component
%     newtheta = theta;
%     newtheta.D{end+1} = theta.D{idx};
%     % initialize the weights
%     newtheta.p(idx) = theta.p(idx) / 2;
%     newtheta.p(end+1) = newtheta.p(idx);
%     % initialize the means
%     Didx = D.component(idx);
%     pp = Didx.llvec(theta.D{idx}, data);
%     [~, I] = sort(pp, 'descend');
%     I = I(1:10);
%     [~, C] = kmeans(data(:,I).', 2);
%     newtheta.D{idx}.mu = C(1,:).';
%     newtheta.D{end}.mu = C(2,:).';
%     % initialize the covariance matrices
%     A = theta.D{idx}.sigma;
%     d = size(A,1);
%     newtheta.D{idx}.sigma = det(A)^(1/d) * eye(d);
%     newtheta.D{end}.sigma = newtheta.D{idx}.sigma;
%   end
%

%% Inner options
% Some estimation functions contain inner estimations that might need to be
% customized separately from the main estimation function. The options for
% the inner estimations can be set in the special field |options.inner|.
% This field can contain usual options as its sub-fields (e.g.
% |options.inner.maxIter|), but some special options like |prevInfo| might
% be disregarded.
%
% Some functions might be more specific about their inner estimation
% options and let you set specific options for each type of inner
% estimation they use. This capability is provided through one more level
% of sub-fields in |options.inner|. For example you can set the maximum
% iterations for the partial inner estimations of the split-and-merge
% algorithms through |options.inner.partial.maxIter|. This should be
% documented in the API reference for the respective function.
%
