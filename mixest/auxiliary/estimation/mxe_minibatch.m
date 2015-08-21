%% |mxe_minibatch|
% *Note:* This is a private function.
%
% Returns a structure used for mini batch optimization
%
% *Syntax*
%
%   [mb, idxMBTrain] = mxe_minibatch(idxAllTrain, options)
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

function [mb, idxMBTrain] = mxe_minibatch(idxAllTrain, options)

    minibatchsize = options.minibatch.size;
    minibatchoverlap = options.minibatch.overlap;
    
    nAllTrain = numel(idxAllTrain);
    idxMiniBatch = 1:min(minibatchsize, nAllTrain);
    idxMBTrain = idxAllTrain(idxMiniBatch);

    mb.next = @next;
    function idxMBTrain = next()
        
        if nAllTrain <= minibatchsize
            idxMiniBatch = 1:nAllTrain;
        else
        
            idxMiniBatchEnd = mod(idxMiniBatch(end) - minibatchoverlap - 1, nAllTrain) + 1;
            if idxMiniBatchEnd >= nAllTrain
                idxMiniBatch = 1:minibatchsize;
            else
                nAvailable = nAllTrain - idxMiniBatchEnd;
                if nAvailable >= minibatchsize
                    idxMiniBatch = idxMiniBatchEnd + 1 : idxMiniBatchEnd + minibatchsize;
                else
                    idxMiniBatch = [idxMiniBatchEnd + 1 : idxMiniBatchEnd + nAvailable, ...
                                    1 : minibatchsize - nAvailable];
                end
            end
            
        end
        idxMBTrain = idxAllTrain(idxMiniBatch);
    end
end