%% |ag_estimatedefault|
% *Note:* This is a private function.
%
% Default estimation function for eg distribution
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Contributors:
%  Reshad Hosseini
%  Yasaman Etesam
%  Pourya H. Zadeh
%
% Change log:
%

function [theta, D, info, options] = ag_estimatedefault(D, data , options)


data = mxe_readdata(data);
weight = data.weight;
N = data.size;
data = data.data;

if nargin < 3
    options = mxe_options();
else
    options = mxe_options(options);
end

if isempty(weight)
    weight = ones(1, N);
end
n = sum(weight);

% digging our parameters from parameter structure
q = size(data, 1);

% computing func for log-likelihood evaluation
func = @(x)(-q/2).*log(x);

% Initializing Optimization
theta0 = options.theta0;
if isempty(theta0)
    M = eye(q);
else
    theta0 = D.init(data);
    M = theta0.sigma;
end

ll_old = -Inf;
if options.penalize
    invLambdaChol = chol(options.penalizertheta.invLambda);
    alpha = options.penalizertheta.alpha;
end


for it = 1:options.maxiter
    R = chol(M);
    logdetM = 2*sum(log(diag(R)));
    Rinv = R \ eye(q); % faster version of inv_triu(R);
    u = weight ./ sum((Rinv.' * data).^2, 1);
  

    ut = sqrt(u);
    dataw = bsxfun(@times, ut, data);
    
    
    M = (dataw * dataw.');
    
    if options.penalize
        Mmid = invLambdaChol * Rinv;
        denominator = sum(Mmid(:).^2);
        M = alpha/denominator*options.penalizertheta.invLambda + M;
        M = q/(sum(u)+1) * M;
    else
        M = q/sum(u) * M;
    end

    %trace(M)
    M = M /trace(M);
    
    fv = func(u);
    if options.penalize
        ll = -(n+alpha)/(2*n)*logdetM + sum(fv)/n + alpha*func(denominator)/n;
    else
        ll =  (-1/2)*logdetM  +  sum(fv)/n;
    end
    ll_diff = ll - ll_old;
    
    if options.verbosity >= 2
        if options.penalize
            fprintf('Penalized Log likelihood = %g , PLL diff= %g \n', ll, ll_diff);
        else
            fprintf('Log likelihood = %g , LL diff= %g \n', ll, ll_diff);
        end
    end
    if ll_diff <= options.tolcostdiff
        break;
    end
    
    ll_old = ll;
end
 %   ll_diff
theta.sigma = M;
if nargout > 2
    info = [];
end
