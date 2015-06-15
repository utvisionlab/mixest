%% Estimation Statistics Structure

%%
% Every iterative estimation function in MixEst, returns some information
% about its progress at each iteration. This information is output as the
% structure array |info| from the estimation function:
%
%   [theta, D, info, options] = estimation_function(...)
%
% Each element in |info| is a structure with fields corresponding to each
% logged info. For instance |info(1).iter| contains zero (since the
% iterations are numbered starting at zero) and |info(end).iter| contains
% the last iteration number. You may also use the |[info.cost]| syntax to
% get all the logged costs for all iterations as a vector. The following
% example plots the cost by iterations:
%
%   plot([info.iter], [info.cost])
%
% You can find information about the fields present in |info| output for
% each estimation function in its documentation. When a Manopt solver is
% used in the estimation, the |info| returned by the Manopt solver is
% returned directly.
%