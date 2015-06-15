%% |mxe_BIC|
% *Note:* This is a private function.
%
% Calculate BICc
%
% *Syntax*
%
%   bic = mxe_BIC(D, data)
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

function bic = mxe_BIC(D, data)
%
%  bic=icBIC(theta,dat)
%
% The function calculates (negative)BIC information criteria for the model 
%
%  Inputs:
%    theta  : theta is mcgsm distribution class
%    Data: A NISDET class variable for Data
%
%  Outputs:
%    aic: the total criteria value

    data = mxe_readdata(data);
    
    n = D.dim();
    N = data.size;
    bic = -n/2*log(N);
end        
