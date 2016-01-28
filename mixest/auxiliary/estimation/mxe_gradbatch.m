%% |mxe_costgrad|
% *Note:* This is a private function.
%
% Returns the Riemannian gradient of cost (negative log-likelihood)
% with respect to distribution parameters calculated on data, considering
% penalization, etc. Used in gradient-based estimation functions.
% 
% Since it is for SGD no databatchsize is necessary
% 
%
% *Syntax*
%
%   grad = mxe_gradbatch(D, theta, data, options)
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

function grad = mxe_gradbatch(D, theta, data, batch_index, options)
% Note: options are NOT optional here. Following fields are required:
%
% * penalize
% * penalizerTheta  (must contain valid penalizer parameters if penalize is true)
% * dataPatchSize
%

    if nargin < 5
        options.sgd.batchnum = 1;
        batch_index = 1;
    end
    % We don't use Manopt stores because it reduces performance due to
    % calculating the gradient in line-search.
    data = mxe_readdata(data, false);
    idx = data.index;
    datamat = data.data;
    data_size = length(idx);
        
    batchnum = options.sgd.batchnum;
    batch_size = floor(data_size / batchnum);
    index_begin = (batch_index-1) * batch_size + 1;
    index_end = min(batch_index*batch_size, data_size);
    if batch_index == batchnum
        index_end = data_size;
    end

    if options.penalize
        % Calculating the penalizer gradient
        egradPen = D.penalizergrad(theta, options.penalizertheta);
    end
    
    egrad = D.llgrad(theta, datamat(:, idx(index_begin:index_end)));
    egrad = D.scaleparam(-1, egrad);
    if options.penalize
        egradPen = D.scaleparam(-1/batchnum, egradPen);
        egrad = D.sumparam(egrad, egradPen);
    end
    
    grad = D.M.egrad2rgrad(theta, egrad);
    
end
