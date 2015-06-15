%% |check_kl|
% *Note:* This is a private function.
%
% Checking kl-divergence derivation
%
% *Syntax*
%
%   check_kl(D, p1, p2)
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

function check_kl(D, p1, p2)

if nargin < 3
    p2 = D.randparam();
end
if nargin < 2
    p1 = D.randparam();
end
nrd = 200000;
disp(['Numerical Checking using ' num2str(nrd) ' datapoints']);
data = D.sample(p1, nrd);

% First checking the entropy derivation
disp(' ');
disp('+++ Checking Entropy -> entropy(theta1) +++');
disp(['       Numerical Entropy:' num2str(-D.ll(p1,data)/nrd)]);
disp(['       Theoretical Entropy:' num2str(D.entropy(p1))]);

% Now checking kl-divergence derivation
disp(' ');
disp('+++ Checking KL-divergence -> kl(theta1||theta2) +++');
disp(['       Numerical KL:' num2str((D.ll(p1,data)-D.ll(p2,data))/nrd)]);
disp(['       Theoretical KL:' num2str(D.kl(p1,p2))]);