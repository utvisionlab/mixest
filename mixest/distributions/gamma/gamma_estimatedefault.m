%% |gamma_estimatedefault|
% *Note:* This is a private function.
%
% Default estimation function for the gamma distribution
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

function [theta, D, info, options] = gamma_estimatedefault(D, data, options)

data = mxe_readdata(data);
weight = data.weight;
datamat = data.data;

if nargin < 3
    options = mxe_options();
else
    options = mxe_options(options);
end

if isempty(options.theta0)
    theta = D.randparam();
else
    theta = options.theta0;
end
theta = D.fullparam(theta);

% Below is the parameters needed for the estimation (sufficient statistics)
if isempty(weight)
    n = data.size;
    par1 = sum( datamat, 2 ) / n;
    par2 = sum(log(datamat), 2) /n;
else
    n = sum(weight);
    par1 = sum(weight .*     datamat , 2 ) / n;
    par2 = sum(weight .* log(datamat) , 2) /n;
end

aa = theta.a;

iterin = 0;
while 1
    iterin = iterin + 1;
    if iterin > 100
        error('Difficulties with estimating the gamma distribution');
    end
    if aa < eps
        aa = eps;
    end
    anew = (1/aa + (par2 - log(par1) + log(aa) - psi(0,aa)) / ...
        (aa * aa * (1/aa - psi(1,aa)))) ^ -1;
    err = abs(anew - aa);
    if err < options.minstepsize
        aa = anew;
        break
    end
    aa = anew;
end
theta.a = aa;
theta.b = par1 / aa;

if nargout > 2
    info = [];
end
