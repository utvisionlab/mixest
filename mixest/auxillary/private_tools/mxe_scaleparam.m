%% |mxe_scaleparam|
% *Note:* This is a private function.
%
% Scaling the parameters
%
% *Syntax*
%
%   param = mxe_scaleparams(scale, param1)
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

function param = mxe_scaleparam(scale, param1)
param = param1;
if isobject(param1)
    param = param1.sumparam(scale, param1);
elseif iscell(param1)
    for k=1:length(param1)
        param{k} = mxe_sumparam(scale, param1{k});
    end
elseif isstruct(param1)
    st=fieldnames(param1);
    for k=1:length(st)
        param.(st{k}) = mxe_sumparam(scale, param1.(st{k}));
    end
elseif isnumeric(param1)
    param = scale * param1;
else ischar(obj)
    param = param1;
end