%% Distribution Parameters Structure
% Since by design, we keep the distribution parameters apart from the
% distribution structure, many functions of the distributions require a
% parameter structure |theta| to be passed in, or output this parameter
% structure. In what follows we describe this structure in detail (well,
% there is not that much detail, actually).
%

%% Structure Description
% |theta| is a structure with field names corresponding to parameter names
% defined in the respective distribution, and values equal to the desired
% parameter values.
%
% *Example*
%
% The parameters for a multi-variate normal distribution are |mu| for the
% mean, and |sigma| for the covariance matrix. So we can construct a
% parameter structure for this distribution like this:
%
%   theta = struct('mu', [1; 2], 'sigma', [2 1; 3 4]);
%
% or equivalently like this:
%
%   theta.mu = [1; 2];
%   theta.sigma = [2 1; 3 4];
%
% Then we can pass |theta| to the distribution functions. For example:
%
%   D = mvnfactory(2);
%   h = D.entropy(theta)
%
%   h =
%       3.8108
%

%% Relation with Manopt Manifolds
% MixEst distribution parameter structures (|theta|) are actually points on
% some specific manifold representing the compound geometry of valid values
% for all of the parameters of the distribution. The manifold is created
% using the |productmanifold| function from the <http://www.manopt.org
% Manopt toolbox> to build a product of the manifolds corresponding to each
% parameter, and therefore you can work with the value of each parameter
% using a field in |theta|.
%
% The Manopt manifold of parameters for a distribution structure |D| can be
% accessed using |D.M|.
%
