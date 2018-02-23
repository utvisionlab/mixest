%% Distribution Structure Common Members
% MixEst distribution structures expose functionality through structure
% fields which are (mostly) function handles. Here is the list of common
% fields available in all distributions (in no special order).

%% |name|
% Distribution name
%
% *Syntax*
%
%   str = D.name()
%
% *Description*
%
% |str = D.name()| returns the name of the distribution |D| as a string.
%
% *Note:* You need to include the parentheses |()| for this to work.
%
% *Example*
%
%   D = mvnfactory(1);
%   D.name()
% 
%  ans =
%       mvn
%

%% |M|
% Distribution parameter manifold.
%
% *Syntax*
%
%   D.M
%
% *Description*
%
% |D.M| contains the <http://manopt.org Manopt> manifold structure
% corresponding to the parameters of the distribution |D|.
%
% *Example*
%
%   D = mvnfactory(2);
%   D.M.name()
% 
%  ans =
%       Product manifold: [mu: Euclidean space R^(2x1)] x [sigma: SPD manifold (2, 2)]
%

%% |dim|
% Parameter-space dimensions
%
% *Syntax*
%
%   dim = D.dim()
%
% *Description*
%
% |dim = D.dim()| returns the dimensions of the parameter space of the
% distribution |D|.
%
% *Note:* You need to include the parentheses |()| for this to work.
%
% *Example*
%
% The parameter space dimensions of a multi-variate normal distribution
% defined over a 2-dimensional data space is six: two for the mean vector
% and four for the 2-by-2 covariance matrix.
%
%   D = mvnfactory(2);
%   dim = D.dim()
% 
%   dim =
%        6
%

%% |datadim|
% Data-space (sample) dimensions
%
% *Syntax*
%
%   n = D.datadim()
%
% *Description*
%
% |n = D.datadim()| returns the dimensions of the data space where the
% distribution |D| is defined.
%
% *Note:* You need to include the parentheses |()| for this to work.
%
% *Example*
%
% The following example creates a multi-variate normal distribution defined
% over a 2-dimensional data space and calls |datadim| to get its data-space
% dimensions.
%
%   D = mvnfactory(2);
%   n = D.datadim()
% 
%   n =
%        2
%

%% |ll|
% Log-likelihood
%
% *Syntax*
%
%   ll = D.ll(theta, data)
%   [ll, store] = D.ll(theta, data, store)
%
% *Description*
%
% |ll = D.ll(theta, data)| returns the total log-likelihood of the
% distribution |D| calculated over |data|, given the distribution
% parameters |theta|.
%
% |[ll, store] = D.ll(theta, data, store)| can be used for caching
% purposes as described <../caching.html here>.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *Example*
%
%   % create a distribution and random parameters
%   D = mvnfactory(1);
%   theta = D.randparam();
%   % generate 1000 random data points
%   data = randn(1, 1000);
%   % find the log-likelihood
%   ll = D.ll(theta, data)
%

%% |llvec|
% Log-likelihood for each data point
%
% *Syntax*
%
%   llvec = D.llvec(theta, data)
%   [llvec, store] = D.llvec(theta, data, store)
%
% *Description*
%
% |llvec = D.llvec(theta, data)| returns the log-likelihood for each 
% point in |data| as a |1-by-N| vector where |N| is the number of data
% points. |theta| is the distribution parameters.
%
% |[llvec, store] = D.llvec(theta, data, store)| can be used for caching
% purposes as described <../caching.html here>.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *Example*
%
%   % create a distribution and random parameters
%   D = mvnfactory(1);
%   theta = D.randparam();
%   % generate 1000 random data points
%   data = randn(1, 1000);
%   % find the log-likelihood vector
%   llvec = D.llvec(theta, data)
%

%% |llgrad|
% Gradient of the log-likelihood, with respect to parameters
%
% *Syntax*
%
%   dll = D.llgrad(theta, data)
%   [dll, store] = D.llgrad(theta, data, store)
%
% *Description*
%
% |dll = D.llgrad(theta, data)| returns the (Euclidean) gradient of the
% log-likelihood (of the distribution |D| calculated over |data|, given the
% distribution parameters |theta|), with respect to parameters, in the
% parameter structure |dll|.
%
% |[dll, store] = D.llgrad(theta, data, store)| can be used for caching
% purposes as described <../caching.html here>.
%
% For information about the parameter input |theta| and output |dll|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *Example*
%
%   % create a distribution and random parameters
%   D = mvnfactory(1);
%   theta = D.randparam();
%   % generate 1000 random data points
%   data = randn(1, 1000);
%   % This is how you should call the functions if you 
%   % need the log-likelihood along with its gradient, in
%   % order to share the common intermediate variables:
%   [ll, store] = D.ll(theta, data);
%   dll = D.llgrad(theta, data, store);
%

%% |llgraddata|
% Gradient of the log-likelihood, with respect to data
%
% *Syntax*
%
%   dld = D.llgraddata(theta, data)
%   [dld, store] = D.llgraddata(theta, data, store)
%
% *Description*
%
% |dld = D.llgraddata(theta, data)| returns the (Euclidean) gradient of the
% log-likelihood (of the distribution |D| calculated over |data|, given the
% distribution parameters |theta|), with respect to each data point, in the
% |n-by-N| vector |dld|, where |n| is the data dimensions and |N| is the
% number of data points.
%
% |[dld, store] = D.llgraddata(theta, data, store)| can be used for caching
% purposes as described <../caching.html here>.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *Example*
%
%   % create a distribution and random parameters
%   D = mvnfactory(1);
%   theta = D.randparam();
%   % generate 1000 random data points
%   data = randn(1, 1000);
%   % This is how you should call the functions if you 
%   % need the log-likelihood along with its gradient, in
%   % order to share the common intermediate variables:
%   [ll, store] = D.ll(theta, data);
%   dld = D.llgraddata(theta, data, store);
%

%% |cdf|
% Cumulative distribution function
%
% *Syntax*
%
%   y = D.cdf(theta, data)
%
% *Description*
%
% |y = D.cdf(theta, data)| returns the cumulative distribution function
% (CDF) of the distribution |D| evaluated at points in |data|, given the
% distribution parameters |theta|. |y| is a |1-by-N| vector where |N| is
% the number of data points.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *Example*
%
%   % create a distribution and random parameters
%   D = mvnfactory(1);
%   theta = D.randparam();
%   % generate 1000 random data points
%   data = randn(1, 1000);
%   % find the CDF on data
%   y = D.cdf(theta, data)
%

%% |pdf|
% Probability density function 
%
% *Syntax*
%
%   y = D.pdf(theta, data)
%
% *Description*
%
% |y = D.pdf(theta, data)| returns the probability density function
% (PDF) of the distribution |D| evaluated at points in |data|, given the
% distribution parameters |theta|. |y| is a |1-by-N| vector where |N| is
% the number of data points.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *Example*
%
%   % create a distribution and random parameters
%   D = mvnfactory(1);
%   theta = D.randparam();
%   % generate 1000 random data points
%   data = randn(1, 1000);
%   % find the PDF on data
%   y = D.pdf(theta, data)
%

%% |sample|
% Generate random samples
%
% *Syntax*
%
%   data = D.sample(theta)
%   data = D.sample(theta, num)
%
% *Description*
%
% |data = D.sample(theta)| generates a random sample drawn from the
% distribution |D|, given the distribution parameters |theta|. |data| is an
% |n-by-1| vector where |n| is the dimensions of the data space.
%
% |data = D.sample(theta, num)| generates |num| random samples. |data| is
% an |n-by-num| matrix.
%
% *Example*
%
% Generate 1000 random samples from an N(0,1) distribution:
% 
%   D = mvnfactory(1);
%   theta = struct('mu', 0, 'sigma', 1);
%   data = D.sample(theta, 1000)
%

%% |randparam|
% Generate random parameters for the distribution
%
% *Syntax*
%
%   theta = D.randparam()
%
% *Description*
%
% |theta = D.randparam()| generates a valid random parameter structure for
% the distribution |D|.
%
% *Note:* You need to include the parentheses |()| for this to work.
%
% For more information about the output |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>.
%
% *Example*
%
% Generate a random parameter structure for a Gaussian distribution:
% 
%   D = mvnfactory(1);
%   theta = D.randparam()
%

%% |init|
% Generate initial parameters appropriate for the given data
%
% *Syntax*
%
%   theta = D.init(data)
%   theta = D.init(data, ...)
%
% *Description*
%
% |theta = D.init(data)| generates suitable parameters to be used as the
% initial point for later estimation on |data|.
%
% |theta = D.init(data, ...)| Some distributions may accept additional
% arguments for the |init| function. This should be mentioned in their
% specific documentation.
%
% For information about the output |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *Example*
%
%   % create a mixture of two Gaussian distributions
%   D = mixturefactory(mvnfactory(1), 2);
%   % generate 1000 random data points
%   data = [randn(1,500), randn(1,500)+5];
%   % find a suitable initialization point
%   options.theta0 = D.init(data);
%   % perform estimation
%   theta = D.estimate(data, options)
%

%% |estimate|
% Estimate distribution parameters to fit data
%
% *Syntax*
%
%   theta = D.estimate(data)
%   theta = D.estimate(data, options)
%   [theta, D] = D.estimate(...)
%   [theta, D, info] = D.estimate(...)
%   [theta, D, info, options] = D.estimate(...)
%
% *Description*
%
% |theta = D.estimate(data)| returns estimated parameters for the
% distribution |D|, using |data|.
%
% |theta = D.estimate(data, options)| utilizes applicable options
% from the |options| structure in the estimation procedure.
%
% |[theta, D] = D.estimate(...)| also returns |D|, the distribution
% structure for which |theta| is applicable. (This is the same as the
% distribution structure |D| from which you called |estimate|, and so it
% should not normally be used. The purpose of including it in the output is
% to maintain compatibility with other estimation functions).
%
% |[theta, D, info] = D.estimate(...)| also returns |info|, a
% structure array containing information about successive iterations
% performed by iterative estimation functions.
%
% |[theta, D, info, options] = D.estimate(...)| also returns the
% effective |options| used, so you can see what default values the function
% used on top of the options you possibly specified.
%
% For information about the output |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>. The
% input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>. You may also want to read about
% <../estimation_options.html |options|> or
% <../estimation_statistics_structure.html |info|> arguments.
%
% *Available Options*
%
% This function supports all the options described in
% <../estimation_options.html estimation options>.
%
% *Returned |info| fields*
%
% The fields present in the returned |info| structure array, depend on the
% solver used (|options.solver|). When a Manopt solver is specified, the
% |info| returned by the Manopt solver is returned directly. For the 'default'
% solver see the documentation of the 'estimatedefault' function for the
% specific distribution. You can read more at our documentation on
% <../estimation_statistics_structure.html estimation statistics
% structure>.
%
% *Example*
%
%   % create a Gaussian distribution
%   D = mvnfactory(1);
%   % generate 1000 random data points
%   data = randn(1,1000) .* 2 + 1;
%   % set some options
%   options.solver = 'cg';
%   options.verbosity = 2;
%   % fit distribution parameters to data
%   theta = D.estimate(data, options)
%

%% |penalizerparam|
% Generate parameter structure for the default penalization function
%
% *Syntax*
%
%   penalizer_theta = D.penalizerparam(data)
%
% *Description*
%
% |penalizer_theta = D.penalizerparam(data)| returns the parameter
% structure for the default penalization function related to the
% distribution |D|, appropriate for |data|.
%
% For more information on special penalization functions and their
% parameters for each distribution, refer to the documentation of that
% distribution.
%
% The input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%

%% |penalizercost|
% Cost penalty of the default penalization function
%
% *Syntax*
%
%   costP = D.penalizercost(theta, penalizer_theta)
%   [costP, store] = D.penalizercost(theta, penalizer_theta, store)
%
% *Description*
%
% |costP = D.penalizercost(theta, penalizer_theta)| returns the cost
% penalty calculated by the default penalization function related to the
% distribution |D|. |theta| is the distribution parameter structure and
% |penalizer_theta| represents the parameter structure for the penalization
% function.
%
% |[costP, store] = D.penalizercost(theta, penalizer_theta, store)| can be
% used for caching purposes as described <../caching.html here>.
%
% The output of this function is used as a regularizer for the cost
% function during parameter estimation when penalization is turned on.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>.
% |penalizer_theta| can be obtained by calling <#15 |penalizerparam|>.
%

%% |penalizergrad|
% Gradient penalty of the default penalization function
%
% *Syntax*
%
%   gradP = D.penalizergrad(theta, penalizer_theta)
%   [gradP, store] = D.penalizergrad(theta, penalizer_theta, store)
%
% *Description*
%
% |gradP = D.penalizergrad(theta, penalizer_theta)| returns the (Euclidean)
% gradient penalty calculated by the default penalization function related
% to the distribution |D|. |theta| is the distribution parameter structure
% and |penalizer_theta| represents the parameter structure for the
% penalization function.
%
% |[gradP, store] = D.penalizergrad(theta, penalizer_theta, store)| can be
% used for caching purposes as described <../caching.html here>.
%
% The output of this function is used as a regularizer for the
% cost-gradient function during parameter estimation when penalization is
% turned on.
%
% For information about the parameter input |theta| and output |gradP|, see
% <../distribution_parameters.html Distribution Parameters Structure>.
% |penalizer_theta| can be obtained by calling <#15 |penalizerparam|>.
%

%% |sumparam|
% Sum of two parameter structures
%
% *Syntax*
%
%   theta = D.sumparam(theta1, theta2)
%
% *Description*
%
% |theta = D.sumparam(theta1, theta2)| calculates the element-by-element
% sum of two parameter structures |theta1| and |theta2| of the distribution
% |D|. |theta1| and |theta2| can also be Euclidean gradients with respect
% to parameter distributions.
%
% For information about the parameter structures |theta|, |theta1| and
% |theta2|, see <../distribution_parameters.html Distribution Parameters
% Structure>.
%
% *Example*
%
% Following example shows how to add two Euclidean gradients of the
% log-likelihood with respect to parameters:
%
%   % create a distribution and random parameter values
%   D = mvnfactory(1);
%   theta1 = D.randparam();
%   theta2 = D.randparam();
%   % generate 1000 random data points
%   data = randn(1, 1000);
%   % find Euclidean gradients of ll and penalizer
%   dll1 = D.llgrad(theta1, data);
%   dll2 = D.llgrad(theta2, data);
%   % sum the gradients
%   dllsum = D.sumparam(dll1, dll2)
%

%% |scaleparam|
% Multiply parameter structure by a scalar
%
% *Syntax*
%
%   theta = D.scaleparam(scalar, theta)
%
% *Description*
%
% |theta = D.scaleparam(scalar, theta)| calculates the product of the
% scalar |scalar| by the parameter structure |theta| of the distribution
% |D|.
%
% For information about the parameter structure |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>.
%
% *Example*
%
% Following example shows how to negate the Euclidean gradient of the
% log-likelihood with respect to parameters:
%
%   % create a distribution and random parameter values
%   D = mvnfactory(1);
%   theta = D.randparam();
%   % generate 1000 random data points
%   data = randn(1, 1000);
%   % find Euclidean gradient of ll
%   dll = D.llgrad(theta, data);
%   % negate the gradient
%   grad = D.scaleparam(-1, dll)
%

%% |sumgrad|
% Sum of two Riemannian gradients
%
% *Syntax*
%
%   rgrad = D.sumgrad(rgrad1, rgrad2, theta)
%
% *Description*
%
% |rgrad = D.sumgrad(rgrad1, rgrad2, theta)| calculates the sum of two
% Riemannian gradients |rgrad1| and |rgrad2| on the tangent space at the
% point |theta| on the parameter manifold of the distribution |D|. |rgrad1|
% and |rgrad2| should be obtained by |D.M.egrad2rgrad| from the
% corresponding Euclidean gradients.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>.
%
% *Example*
%
% Following example shows how to add the Riemannian gradients of the
% penalizer and log-likelihood at some parameter value:
%
%   % create a distribution and random parameter values
%   D = mvnfactory(1);
%   theta = D.randparam();
%   % generate 1000 random data points
%   data = randn(1, 1000);
%   % find Euclidean gradients of ll and penalizer
%   dll = D.llgrad(theta, data);
%   penalizer_theta = D.penalizerparam(data);
%   gradP = D.penalizergrad(theta, penalizer_theta);
%   % convert the gradients to Riemannian
%   rdll = D.M.egrad2rgrad(theta, dll);
%   rgradP = D.M.egrad2rgrad(theta, gradP);
%   % sum the gradients
%   grad = D.sumgrad(rdll, rgradP, theta)
%

%% |scalegrad|
% Multiply Riemannian gradient by a scalar
%
% *Syntax*
%
%   rgrad = D.scalegrad(scalar, rgrad, theta)
%
% *Description*
%
% |rgrad = D.scalegrad(scalar, rgrad, theta)| calculates the product of
% the scalar |scalar| by the Riemannian gradient |rgrad|, on the tangent
% space at the point |theta| on the parameter manifold of the distribution
% |D|. |rgrad| should be obtained by |D.M.egrad2rgrad| from the
% corresponding Euclidean gradient.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>.
%
% *Example*
%
% Following example shows how to negate the Riemannian gradient of the
% log-likelihood at some parameter value:
%
%   % create a distribution and random parameter values
%   D = mvnfactory(1);
%   theta = D.randparam();
%   % generate 1000 random data points
%   data = randn(1, 1000);
%   % find Euclidean gradient of ll
%   dll = D.llgrad(theta, data);
%   % convert the gradient to Riemannian
%   rdll = D.M.egrad2rgrad(theta, dll);
%   % negate the gradient
%   grad = D.scalegrad(-1, rdll, theta)
%

%% |entropy|
% Calculate entropy
%
% *Syntax*
%
%   h = D.entropy(theta)
%
% *Description*
%
% |h = D.entropy(theta)| calculates the entropy of the distribution |D|
% given parameters |theta|.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>.
%
% *Example*
%
%   % create a distribution and random parameter values
%   D = mvnfactory(2);
%   theta = D.randparam();
%   % find the entropy
%   h = D.entropy(theta)
%

%% |kl|
% Calculate Kullback–Leibler divergence
%
% *Syntax*
%
%   kl = D.kl(theta1, theta2)
%
% *Description*
%
% |kl = D.kl(theta1, theta2)| calculates the KL divergence between the
% distribution |D| given parameters |theta1| and  the same distribution
% given parameters |theta2|.
%
% For information about the parameter inputs |theta1| and |theta2|, see
% <../distribution_parameters.html Distribution Parameters Structure>.
%
% *Example*
%
%   % create a distribution and two random parameter structures
%   D = mvnfactory(2);
%   theta1 = D.randparam();
%   theta2 = D.randparam();
%   % find the KL divergence
%   kl = D.kl(theta1, theta2)
%

%% |AICc|
% Calculate corrected Akaike information criterion without the likelihood
% term
%
% *Syntax*
%
%   aicc = D.AICc(data)
%
% *Description*
%
% |aicc = D.AICc(data)| calculates the corrected Akaike information
% criterion (AICc) of the distribution |D| for the given |data|, without
% the likelihood term. This should be added to the <#5 log-likelihood> to
% obtain the full AICc (multiplied by -1/2).
%
% The input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *Example*
%
%   % create a distribution and random parameters
%   D = mvnfactory(2);
%   theta = D.randparam();
%   % generate 1000 random data points
%   data = randn(2, 1000);
%   % find the AICc
%   aicc = D.ll(theta, data) + D.AICc(data)
%

%% |BIC|
% Calculate Bayesian information criterion without the likelihood term
%
% *Syntax*
%
%   bic = D.BIC(data)
%
% *Description*
%
% |bic = D.BIC(data)| calculates the Bayesian information criterion (BIC)
% of the distribution |D| for the given |data|, without the likelihood
% term. This should be added to the <#5 log-likelihood> to obtain the full
% BIC (multiplied by -1/2).
%
% The input argument |data| is described in <../data_input.html Data Input
% Argument to Functions>.
%
% *Example*
%
%   % create a distribution and random parameters
%   D = mvnfactory(2);
%   theta = D.randparam();
%   % generate 1000 random data points
%   data = randn(2, 1000);
%   % find the BIC
%   bic = D.ll(theta, data) + D.BIC(data)
%

%% |display|
% Display parameter values
%
% *Syntax*
%
%   str = D.display(theta)
%   D.display(theta)
%
% *Description*
%
% |str = D.display(theta)| returns a string containing information about
% the parameter values in |theta|.
%
% |D.display(theta)| displays the distribution name along with information
% about the parameter values in |theta|.
%
% For information about the parameter input |theta|, see
% <../distribution_parameters.html Distribution Parameters Structure>.
%
% *Example*
%
%   % create a distribution and parameter structure
%   D = mvnfactory(1);
%   theta = struct('mu', [0;0], 'sigma', [2 1;3 4]);
%   % display parameter values
%   D.display(theta)
%
%  ans =
%       mvn distribution parameters:
%       mean (2-by-1): [0;0]
%       covariance (2-by-2): [2 1;3 4]
%

%% |selfsplit|
% Calculate the initial parameters for two splitted distributions to be
% substituted for current distribution (used in split-and-merge mixture
% estimation)
%
% *Syntax*
%
%   [value1, value2] = D.selfsplit(theta, param_name)
%   [value1, value2] = D.selfsplit(theta, param_name, method)
%   [value1, value2] = D.selfsplit(theta, param_name, method, data)
%   [value1, value2, store] = D.selfsplit(theta, param_name, method, data, store)
%   [value1, value2, store, mixture_store] = D.selfsplit(theta, param_name, method, data, store, mixture_D, mixture_theta, idx, mixture_store)
%

%% |selfmerge|
% Calculate the initial parameters for a merged distribution to be
% substituted for two distributions (used in split-and-merge mixture
% estimation)
%
% *Syntax*
%
%   value = D.selfmerge(theta1, theta2, param_name, w1, w2)
%   value = D.selfmerge(theta1, theta2, param_name, w1, w2, method)
%   value = D.selfmerge(theta1, theta2, param_name, w1, w2, method, data)
%   [value, store] = D.selfmerge(theta1, theta2, param_name, w1, w2, method, data, store)
%   [value, store, mixture_store] = mvn_selfmerge(theta1, theta2, param_name, w1, w2, method, data, store, mixture_D, mixture_theta, idx1, idx2, mixture_store)
%
