%% |mxe_readweight|
% *Note:* This is a private function.
%
% Read only the data weights from the given data argument
%
% *Syntax*
%
%   weight = mxe_readweight(data)
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

function weight = mxe_readweight(data)
    
    weight = [];
    if isstruct(data)
        if isfield(data, 'weight')
            weight = data.weight;
        end
    end
end
