%% |objlen|
% *Note:* This is a private function.
%
% Calculate the number of numerical elements in obj
%
% *Syntax*
%
%   l = objlen(obj)
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

function l = objlen(obj)

if isobject(obj)
    l = obj.objlen();
elseif iscell(obj)
    l=0;
    for k=1:length(obj)
        l=l+objlen(obj{k});
    end
elseif isstruct(obj)
    st=fieldnames(obj);
    l=0;
    for k=1:length(st)
        l=l+objlen(obj.(st{k}));
    end
elseif isnumeric(obj)
    l=numel(obj);
elseif ischar(obj)
    l=0;
end