%% |vmf_estimatedefault|
% *Note:* This is a private function.
%
% Default estimation function for the von Mises-Fisher distribution
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

function [theta, D, info, options] = vmf_estimatedefault(D, data, options)

data = mxe_readdata(data);
weight = data.weight;
N = data.size;
dim = data.dim;
data = data.data;

if isempty(weight)
    weight = ones(1, N);
end

n = sum(weight);

mu     = sum(bsxfun(@times, weight, data), 2);
norm_mu = norm(mu);
theta.mu = mu / norm_mu;

rbar  = norm_mu/n;

theta.kappa = (rbar*dim - rbar^3)/(1-rbar^2);

if nargout > 2
    info = [];
end
