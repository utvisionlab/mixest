%% |mxe_crossvalidation|
% *Note:* This is a private function.
%
% Construct a structure used for cross-validation
%
% *Syntax*
%
%   [cv, idxTrain] = mxe_crossvalidation(data, cvfraction, cvtoliter, costfun)
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

function [cv, idxTrain] = mxe_crossvalidation(data, cvfraction, cvtoliter, costfun)

    data = mxe_readdata(data, false);
    idxAll = data.index;
    nAll = data.size;
    datamat = data.data;

    % split train and cross-validation data
    traindatafraction = 1 - cvfraction;
    cv.nTrain = ceil(traindatafraction .* nAll);
    idxTrain = idxAll(1:cv.nTrain);
    cv.idxCv = idxAll(cv.nTrain+1:nAll);
    dataCv = struct('data',datamat, 'index',cv.idxCv, 'tag','cv');

    % init cv.bestResult
    cv.bestResult.cost = Inf;
    
    

    cv.statsfun = @statsfun;
    function stats = statsfun(theta, stats)
        
        stats.cvCost = costfun(theta, dataCv);
    end

    cv.stopfun = @stopfun;
    function stopnow = stopfun(theta, info, last, first)
    % stop function must be executed after statsfun
        
        if nargin < 4
            first = 1;
        end
        if last < first
            return
        end
        
        stopnow = false;
        stats = info(last);
        iter = stats.iter;
        cvCost = stats.cvCost;
        
        % testCost must be less than at least one of its recent values
        if cvCost < cv.bestResult.cost
            % update cv.bestResult
            cv.bestResult.iter = iter;
            cv.bestResult.theta = theta;
            cv.bestResult.cost = cvCost;
            
        elseif iter > cv.bestResult.iter + cvtoliter
            % no improvement after cvtoliter steps => stop
            stopnow = true;
        end
           
    end

    cv.hasFinalResult = @hasFinalResult;
    function flag = hasFinalResult()
    % always check this function before calling cv.finalResult
    % (maybe stopfun is not called, e.g. considering options.miniter)
    
        flag = isfield(cv.bestResult, 'theta');
    end

    cv.finalResult = @finalResult;
    function [theta, cvCost] = finalResult()
        
        % return the results preceding the cost increase
        theta = cv.bestResult.theta;
        if nargout > 1
            cvCost = cv.bestResult.cost;
        end
    end

end