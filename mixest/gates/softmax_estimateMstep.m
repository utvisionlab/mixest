%% |softmax_estimateMstep|
% *Note:* This is a private function.
%
% Default estimation function for the mixtures of experts. Implements the
% M-step in Expectation Maximization (EM) algorithm.
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

function theta = softmax_estimateMstep(D, data, options)


    datadim = D.datadim();
    % Change of ll and llgrad to use in optimization
    
    %%

    function store = weighting_intermediate_params(theta, data, store)
    % calculate intermediate parameters for ll #1
        if ~isfield(store, 'dataTW')
            data = mxe_readdata(data);
            n = data.size;
            data = data.data;
            store.dataTW = theta.W * [data(1:datadim,:); zeros(1,n)];
        end
    end        


    function store = weighting_intermediate_params2(theta, data, store)
        % calculate intermediate parameters for ll #2
        
        
        store = weighting_intermediate_params(theta, data, store);
        
        if ~isfield(store, 'logpgate')
            logpgate = [zeros(1, size(store.dataTW,2)); store.dataTW];
            store.logpgatesum = logsumexp(logpgate);
        end
            
        if ~isfield(store, 'logpcond')
            data = mxe_readdata(data);
            data = data.data;
            store.logpcond = logpgate .* data(datadim+1:end,:);
            store.logpcondsum = sum(store.logpcond);
        end
    end

%%

    D.llvec = @llvec;
    function [llvec, store] = llvec(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        weight = mxe_readweight(data);
        
        store = weighting_intermediate_params2(theta, data, store);
        
        llvec = store.logpcondsum - store.logpgatesum;
        
        if ~isempty(weight)
            llvec = llvec .* weight;
        end
        
    end

%%

    D.ll = @ll;
    function [ll, store] = ll(theta, data, store)
        
        if nargin < 3
            store = struct;
        end
        
        
        [llvec, store] = D.llvec(theta, data, store);
         
        ll = sum(llvec, 2);

    end

%%

    D.llgrad = @llgrad;
    function [dll, store] = llgrad(theta, data, store)
        
        data = mxe_readdata(data);
        weight = data.weight;
        n = data.size;
        data = data.data;
        
        if nargin < 3
            store = struct;
        end
        
        store = weighting_intermediate_params2(theta, data, store);
        
        pgate = exp(bsxfun(@minus, store.dataTW, store.logpgatesum));
        
        part = data(datadim+2:end,:) - pgate;
        
        datap = [data(1:datadim,:); zeros(1,n)];
        
        if ~isempty(weight)
            datap = bsxfun(@times, weight, datap);
        end

        dll.W = part * datap.';
        
    end

    D = mxe_addsharedfields(D);
    
    theta = D.estimate(data, options);
  
end
