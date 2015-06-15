%% |mxe_AICc|
% *Note:* This is a private function.
%
% Calculate AICc
%
% *Syntax*
%
%   aicc = mxe_AICc(D, data)
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

function aicc = mxe_AICc(D, data)

    data = mxe_readdata(data);
    
    n = D.dim();
    N = data.size;
    aicc = N/(N-n-1)*n;
end
