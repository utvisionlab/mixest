%% |mxe_costgrad|
% *Note:* This is a private function.
%
% Returns the cost (negative log-likelihood) and its Riemannian gradient
% with respect to distribution parameters calculated on data, considering
% penalization, etc. Used in gradient-based estimation functions.
%
% *Syntax*
%
%   [cost, grad] = mxe_costgrad(D, theta, data, options)
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

function [cost, grad] = mxe_costgrad(D, theta, data, options)
% Note: options are NOT optional here. Following fields are required:
%
% * penalize
% * penalizerTheta  (must contain valid penalizer parameters if penalize is true)
% * dataPatchSize
%

    % We don't use Manopt stores because it reduces performance due to
    % calculating the gradient in line-search.
    store = struct;
    %TODO make a convention for a specific field to contain theta-related intermediate parameters
    
    data = mxe_readdata(data, false);
    index = data.index;
    N = data.size;
        
    if options.penalize && ~isempty(options.penalizertheta)
        % Calculating the penalizer cost and gradient
        [costPen, store] = D.penalizercost(theta, options.penalizertheta, store);
        if nargout > 1
            %[egradPen, store] = D.penalizergrad(theta, options.penalizertheta, store);
            egradPen = D.penalizergrad(theta, options.penalizertheta, store);
        end
    end
    
    % load data patch-by-patch and sum the log-likelihood
    if options.datapatchsize < N
        nPatches = floor((N - 1) / options.datapatchsize) + 1;
        %storeCopy = store;
        %if ~isfield(store, 'patches')
        %    store.patches = cell(nPatches,1);
        %end
    else
        nPatches = 1;
        %patchStore = store;
    end
    
    datapatch = data;
    indexb = 0;
    for mini = 1:nPatches
        indexe = min(N, indexb + options.datapatchsize);
        datapatch.index = index(indexb+1:indexe);
        
        %if options.datapatchsize < n
            %patchStore = store.patches{mini};
        %    if ~isstruct(patchStore)
        %        patchStore = storeCopy;
        %    end
        %end
        
        %[lik, patchStore] = D.ll(theta, datapatch, patchStore);
        [lik, patchStore] = D.ll(theta, datapatch);
        if options.regularize
            costR = D.regcost(theta, datapatch, patchStore);
        end
        if nargout > 1
            if options.regularize
                egrad1 = D.reggrad(theta, datapatch, patchStore);
            else
                egrad1 = D.llgrad(theta, datapatch, patchStore);
            end
            %[egrad1, patchStore] = D.llgrad(theta, datapatch, patchStore);
        end

        %if options.datapatchsize < n
        %    store.patches{mini} = patchStore;
        %end
        
        indexb = indexe;
        if mini == 1
            ll = lik;
            if options.regularize
                costReg = costR;
            end
            if nargout > 1
                egrad = egrad1;
            end
        else
            ll = ll + lik;
            if options.regularize
                costReg = costReg + costR;
            end
            if nargout > 1
                egrad = D.sumparam(egrad, egrad1);
            end
        end
    end

    if options.penalize && ~isempty(options.penalizertheta)
        ll = ll + costPen;
        if nargout > 1
            egrad = D.sumparam(egrad, egradPen);
        end
    end
    if options.regularize
        ll = ll + costReg;
    end
    
    cost = (-1/N) * ll;
    if nargout > 1
        grad = D.M.egrad2rgrad(theta, egrad);
        grad = D.M.lincomb(theta, -1/N, grad);
    end
end
