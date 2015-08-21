%% |mxe_sumparam|
% *Note:* This is a private function.
%
% Adding two parameters
%
% *Syntax*
%
%   param = mxe_addparams(param1, param2)
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

function param = mxe_sumparam(param1, param2)
param = param1;
if isobject(param1)
    param = param1.sumparam(param1, param2);
elseif iscell(param1)
    for k=1:length(param1)
        param{k} = mxe_sumparam(param1{k}, param2{k});
    end
elseif isstruct(param1)
    st=fieldnames(param1);
    for k=1:length(st)
        param.(st{k}) = mxe_sumparam(param1.(st{k}), param2.(st{k}));
    end
elseif isnumeric(param1)
    param = param1 + param2;
else ischar(obj)
    param = param1;
end