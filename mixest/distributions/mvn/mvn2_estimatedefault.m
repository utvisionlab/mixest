%% |mvn2_estimatedefault|
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

function [theta, D, info, options] = mvn2_estimatedefault(D, data, options)

data = mxe_readdata(data);
weight = data.weight;
N = data.size;
data = data.data;
data = [data; ones(1,N)];

if isempty(weight)
    weight = ones(1, N);
end
if nargin > 2 && options.penalize
    nu = options.penalizertheta.nu;
    invLambda = options.penalizertheta.invLambda;
    d = size(data,1);
end
n = sum(weight);


theta.sigmat = bsxfun(@times, data, weight) * data' / n;

if nargin > 2 && options.penalize
    theta.sigmat = (invLambda + theta.sigmat *n) / (nu+n+d+1);
end
if nargout > 2
    info = [];
end
