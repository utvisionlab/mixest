%% |ecd_eg_estimatedefault|
% *Note:* This is a private function.
%
% Default estimation function for eg distribution
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

function [theta, D, info, options] = ecd_eg_estimatedefault(D, radialD, data, options)

data = mxe_readdata(data);
weight = data.weight;
N = data.size;
data = data.data;

if nargin < 4
    options = mxe_options();
else
    options = mxe_options(options);
end

if isempty(weight)
    weight = ones(1, N);
end
n = sum(weight);

% digging our parameters from parameter structure
a = radialD.a;
b = radialD.b;
q = size(data, 1);

% computing func for log-likelihood evaluation
func = @(x)2/b*(x)+2*(q/2-a).*log(x);

% Whitenning the data
C = bsxfun(@times, data, weight) * data' / n;
R = chol(C); % C=R' R ( R upper trangular)
Rinv = R \ eye(q); % faster version of inv_triu(R);
data = Rinv' * data; % whitenned data
logdetR = sum(log(diag(R)));

ll_old = -Inf;
for it = 1:options.maxiter
    if a < q/2
        % First optimization approach
        if it == 1
            theta0 = options.theta0;
            if isempty(theta0) 
                M = eye(q);
            else
                M = (a*b)/q * Rinv' * theta0.sigma * Rinv;
            end
        end
        [M, u, logdetM] = updateConcave(a, b, M, data, weight, n);
        %alpha = sum(u)/q/n;
        %M = M * alpha;
    else
        % Second optimization approach
         if it == 1
            theta0 = options.theta0;
            if isempty(theta0) 
                iM = eye(q);
            else
                Rs = chol(theta0.sigma);
                Rsinv = Rs \ eye(q);
                Rsinv = R * Rsinv;
                iM = q/(a*b) * (Rsinv * Rsinv.');
            end
        end
        [iM, u, logdetiM] = updateConvex(a, b, iM, data, weight, n);
        logdetM = -logdetiM;
    end
    
    fv = func(u);
    ll =  logdetM + 2*logdetR +  sum(fv)/n;
    ll_diff = ll - ll_old;
    
    if options.verbosity >= 2
        fprintf('Log likelihood = %g , LL diff= %g \n', ll, ll_diff);
    end
    
    if ll_diff <= options.tolcostdiff && ll > ll_old
        break;
    end
    
    ll_old = ll;
end

if a >= q/2
    Mchol = chol(iM); % M=R' R ( R upper trangular)
    Michol = Mchol \ eye(q);
    M = Michol * Michol.';
end

theta.sigma = R.' * M * R;

if nargout > 2
    info = [];
end

function [Mout, u, logdetM] = updateConcave(a, b, M, data, weight, n)
% The update for concave case a < q/2
% n is sum of weights
% Data is whitened
% Mout = (q-2a)/n \sum w_i [x_i x_i^T] / (x_i^T M^-1 x_i) + 2/b I
q = size(data,1);
%
if ~isequal(M, eye(q))
    R = chol(M);
    logdetM = 2*sum(log(diag(R)));
    Rinv = R \ eye(q); % faster version of inv_triu(R);
    alpha = norm(Rinv,'fro')^2 / ( a * b);
    u = weight ./ sum((Rinv.' * data).^2, 1);
else
    logdetM = 0;
    u = weight ./ sum(data.^2, 1);
    alpha = 1;
end
data = bsxfun(@times, data, sqrt(u));
Mout = (data * data.');
Mout = (q-2*a)/n * alpha * Mout + 2/b * eye(q);

function [iMout, u, logdetiM] = updateConvex(a, b, iM, data, weight, n)
% The update for convex case a > q/2
% n is sum of weights
% Data is whitened
% Mout^-1 = b/2*(q-2a)/n \sum w_i [M^-0.5 x_i x_i^T M^-0.5] / (x_i^T M^-1
% x_i) + b/2 I
q = size(data,1);
%
if ~isequal(iM, eye(q))
    [U,T] = schur(iM,'complex');
    Ts = diag(sqrt(diag(T))); 
    sM = U * Ts * U';
    logdetiM = sum(log(diag(T)));
    u = weight ./ sum((sM * data).^2, 1);
else
    logdetiM = 0;
    u = weight ./ sum(data.^2, 1);
end
data = bsxfun(@times, data, sqrt(u));
iMout = (data * data.');
if ~isequal(iM, eye(q))
    iMout = b*(a-q/2)/n * sM * iMout * sM + b/2 * eye(q);
else
    iMout = b*(a-q/2)/n * iMout + b/2 * eye(q);
end