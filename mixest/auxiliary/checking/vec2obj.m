%% |vec2obj|
% *Note:* This is a private function.
%
% Convert the vector |vec| into an object like the sample object |obj|
%
% *Syntax*
%
%   obj = vec2obj(obj, vec)
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

function obj = vec2obj(obj, vec)

if isobject(obj)
    obj = obj.vec2obj(vec);
elseif iscell(obj)
    beginIt = 0;
    for k = 1:length(obj)
        l = objlen(obj{k});
        if l > 0
            obj{k} = vec2obj( obj{k}, vec(beginIt+1:beginIt+l) );
        end
        beginIt = beginIt + l;
    end
elseif isstruct(obj)
    st = fieldnames(obj);
    beginIt = 0;
    for k = 1:length(st)
        l = objlen(obj.(st{k}));
        if l > 0
            obj.(st{k}) = vec2obj( obj.(st{k}), vec(beginIt+1:beginIt+l) );
        end
        beginIt = beginIt + l;
    end
elseif isnumeric(obj)
    obj = reshape(vec, size(obj));
end