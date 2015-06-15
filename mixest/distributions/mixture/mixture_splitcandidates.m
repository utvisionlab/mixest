%% |mixture_splitcandidates|
% *Note:* This is a private function.
%
% Sort mixture components for splitting, based on given split method
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

function idx = mixture_splitcandidates(D, theta, data, options, n)

    num = D.num();
    
    if nargin<4 || ~isfield(options, 'sm')
        options = mxe_options();
    end
    if ~ischar(options.sm.splitcriterion)
        error('mixture.splitcandidates: Invalid options.splitcriterion')
    end
    if nargin < 5
        n = num;
    end

    switch options.sm.splitcriterion

        case 'kl'
            % Find the one with maximum kl divergence between
            % model distribution and empirical distribution
            h = loc_kl(D, theta, data);
            flagmax = true;

        case 'mll'
            % Find the one with minimum mean local likelihood
            h = mean_loc_ll(D, theta, data);
            flagmax = false;

        case 'entropy'
            % Find the one with maximum entropy
            h = zeros(num, 1);
            for k = 1 : num
                Dk = D.component(k);
                h(k) = Dk.entropy(theta.D{k});
            end
            flagmax = true;

        case 'rand'
            % Random Split
            h = rand(num, 1);

        otherwise
            error('mixture.splitcandidates: Split criterion ''%s'' not recognized', options.sm.splitcriterion)
    end

    % flagmax=true: sort descending (find maximum), flagmax=false: sort ascending (find minimum)
    if flagmax
        [unused, idx] = sort(h, 'descend'); %#ok<ASGLU>
    else
        [unused, idx] = sort(h); %#ok<ASGLU>
    end
    
    if n < num
        idx = idx(1:n);
    end
end


function h = loc_kl(D, theta, data)
% Using local kl divergence between distribution and empirical distribution
% used for splitting

    data = mxe_readdata(data);
    weight = data.weight;
    N = data.size;
    data = data.data;

    if isempty(weight)
        weight = 1;
    end

    % Initialize different variables
    hll = zeros(D.num(), N);
    hllnl = zeros(D.num(), N);

    % Calculate the log-likelihood of data goes to different clusters
    for k = 1:D.num()
        Dk = D.component(k);
        hllnl(k,:) = Dk.llvec(theta.D{k}, data);
        hll(k,:) = log(theta.p(k)) + hllnl(k,:);
    end

    % Based on log-likelihoods calculate the total log-likelihood and weights
    h_sumX = logsumexp(hll,1);

    % Following implements logarithm of above expression
    hXlog = bsxfun(@plus, bsxfun(@minus, hll, h_sumX), log(weight));  

    % Using hX we can compute empirical distribution
    %   Pemp = bsxfun(@rdivide, hX, sum(hX,2));
    % Following two lines implement the logarithm of the above expression
    hXlog_sum = logsumexp(hXlog, 2);
    logPemp = bsxfun(@minus, hXlog, hXlog_sum);
    Pemp = exp(logPemp);

    % Now based on empirical distribution compute minus entropy and cross ent.
    mEnt = sum(Pemp.*logPemp, 2); 

    crossEnt = sum(Pemp.*hllnl, 2);

    h = mEnt - crossEnt;
end


function h = mean_loc_ll(D, theta, data)
% calculating mean local likelihood necessary for splitting

    data = mxe_readdata(data);
    weight = data.weight;
    N = data.size;
    data = data.data;

    if isempty(weight)
        weight = 1;
    end

    % Initialize different variables
    hll = zeros(D.num(), N);

    % Calculate the log-likelihood of data goes to different clusters
    for k = 1:D.num()
        Dk = D.component(k);
        hll(k,:) = log(theta.p(k)) + Dk.ll(theta.D{k}, data, 'vec');
    end

    % Based on log-likelihoods calculate the total log-likelihood and wrights
    h_sumX = logsumexp(hll,1);
    hX = bsxfun(@times, exp(bsxfun(@minus, hll, h_sumX)), weight);   
    hsumCol = sum(hX.*hll,2) ./ sum(hX,2);
    h = hsumCol - log(theta.p);
end
