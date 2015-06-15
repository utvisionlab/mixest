%% |mixture_mergecandidates|
% *Note:* This is a private function.
%
% Sort mixture components for merging, based on given merge method
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

function [idx1, idx2] = mixture_mergecandidates(D, theta, data, options, n)

    num = D.num();
    max_count = num*(num-1)/2;
    
    if nargin<4 || ~isfield(options, 'sm')
        options = mxe_options();
    end
    if ~ischar(options.sm.mergecriterion)
        error('mixture.mergecandidates: Invalid options.mergecriterion')
    end
    if nargin < 5
        n = max_count;
    end
    
    switch options.sm.mergecriterion

        case 'kl'
            % Find using symmetric kl-divergence
            mat = -Inf(num);
            for k1 = 1 : num-1
                for k2 = k1+1 : num
                    D1 = D.component(k1);
                    D2 = D.component(k2);
                    mat(k1,k2) = -(D1.kl(theta.D{k1}, theta.D{k2})+...
                      D2.kl(theta.D{k2}, theta.D{k1})) / 2;  
                end
            end

        case 'overlap'
            % Find using posterior inner product (distribution overlaps)
            mat = post_inner_prod(D, theta, data);

        case 'rand'
            % Random Merge
            mat = triu(rand(num), 1);

        otherwise
            error('mergecandidates: Merge criterion ''%s'' not recognized', options.sm.mergecriterion)
    end
    
    n = min(n, max_count);
    [unused, idx] = sort(mat(:), 'descend'); %#ok<ASGLU>
    idx = idx(1:n); % Note: here, we're removing the zero elements also
    [idx1, idx2] = ind2sub([num,num], idx);
end

function innerprod = post_inner_prod(D, theta, data)
% This function implements posterior inner product

    h = D.weighting(theta, data);

    num = D.num();
    innerprod = zeros(num);
    for k1 = 1 : num-1
        for k2 = k1+1 : num
            innerprod(k1,k2) = h(k1,:) * h(k2,:)'/...
                sqrt((h(k1,:) * h(k1,:)')*(h(k2,:) * h(k2,:)'));
        end
    end
end

