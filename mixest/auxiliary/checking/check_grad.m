%% |check_grad|
% *Note:* This is a private function.
%
% Checking gradient derivation of log-likelihood
%
% *Syntax*
%
%   check_grad(D, theta, data)
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

function check_grad(D, theta, data)

if nargin < 3
    nrd = 10000;
    data = D.sample(theta, nrd);
else
    data = mxe_readdata(data);
    nrd = data.size;
    data = data.data;
end

disp(['Numerical Checking using ' num2str(nrd) ' datapoints']);

% First checking the entropy derivation
disp(' ');
disp('+++ Checking Gradient -> grad(theta) +++');
x0 = obj2vec(theta);
gnum = finidiff(x0, @(x)fun(D, theta, x, data));
[f, gth] = fun(D, theta, x0, data);

disp(['       Numerical Gradient: ' num2str(gnum)]);
disp(['       Theoretical Gradient: ' num2str(gth')]);
disp(['       Maximum Difference: ',num2str(max(abs(gnum-gth')))]);


function g = finidiff(x, fun)

f = fun(x);
len = length(x);
delf = zeros(1,len);
mu = 1e-6*(1+norm(x))/len;
for k = 1:len
    e_k = zeros(len,1);
    e_k(k) = 1;
    delf(k) = fun(x + mu*e_k);
end
g = (delf-f)/mu;


function [y, dy] = fun(D, theta, x, data)

px = vec2obj(theta, x);
if nargout>1
    [y, store] = D.ll(px, data);
    dd = D.llgrad(px, data, store);
    dy = obj2vec(dd);
else
    y = D.ll(px, data);
end