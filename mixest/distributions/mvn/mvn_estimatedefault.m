%% |mvn_estimatedefault|
% *Note:* This is a private function.
%
% Default estimation function for the multi-variate normal distribution
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

function [theta, D, info, options] = mvn_estimatedefault(D, data, options)

data = mxe_readdata(data);
weight = data.weight;
N = data.size;
data = data.data;

if isempty(weight)
    weight = ones(1, N);
end
if nargin > 2 && options.penalize
    kappa = options.penalizertheta.kappa;
    nu = options.penalizertheta.nu;
    mu = options.penalizertheta.mu;
    invLambda = options.penalizertheta.invLambda;
    d = size(data,1);
end
n = sum(weight);

theta.mu = sum(bsxfun(@times, data, weight),2)/n;

data = bsxfun(@minus, data, theta.mu);

theta.sigma = bsxfun(@times, data, weight) * data' / n;

if nargin > 2 && options.penalize
    mat = (kappa*n)/(kappa+n) * ((theta.mu-mu) * (theta.mu-mu).');
    theta.sigma = (invLambda + mat +  theta.sigma *n) / (nu+n+d+1);
    theta.mu = (n * theta.mu + kappa * mu) / (n+kappa);
end
if nargout > 2
    info = [];
end
