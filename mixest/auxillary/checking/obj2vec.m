%% |obj2vec|
% *Note:* This is a private function.
%
% Convert the object into a vector
%
% *Syntax*
%
%   vecOut = obj2vec(obj)
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

function vecOut = obj2vec(obj)

if isobject(obj)
    vecOut = obj.obj2vec();
elseif iscell(obj)
    L=objlen(obj);
    vecOut=zeros(L,1);
    beginIt=0;
    for k=1:length(obj)
        l=objlen(obj{k});
        if l>0
            vecOut( beginIt+1:beginIt+l)=obj2vec(obj{k});
        end
        beginIt=beginIt+l;
    end
elseif isstruct(obj)
    st=fieldnames(obj);
    L=objlen(obj);
    vecOut=zeros(L,1);
    beginIt=0;
    for k=1:length(st)
        l=objlen(obj.(st{k}));
        if l>0
            vecOut(beginIt+1:beginIt+l)=obj2vec( obj.(st{k}) );
        end
        beginIt=beginIt+l;
    end
elseif isnumeric(obj)
    vecOut=obj(:);
else ischar(obj)
    vecOut=[];
end