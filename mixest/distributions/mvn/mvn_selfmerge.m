%% |mvn_selfmerge|
% *Note:* This is a private function.
%
% Calculate the initial parameter values for a merged distribution to be
% substituted for two distributions. Each parameter should be requested by
% a separate call with an optional specific method.
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

function [value, store, mixture_store] = ...
    mvn_selfmerge(D, theta1, theta2, param_name, w1, w2, method, data, store, ...
    mixture_D, mixture_theta, idx1, idx2, mixture_store) %#ok<INUSL>

    if nargin < 6
        method = 'default';
    end
    if nargin < 8
        store = struct;
    end
    
    switch param_name

        case 'mu'
%%
            switch method
                case 'default'
                    if isfield(store, 'selfmerge_theta')
                        theta = store.selfmerge_theta;
                    else
                        theta = D.sumparam(D.scaleparam(w1/(w1+w2),theta1), D.scaleparam(w2/(w1+w2),theta2));
                        store.selfmerge_theta = theta;
                    end
                    value = theta.mu;

                otherwise
                    error('mvn.selfmerge: Method ''%s'' not recognized for parameter ''%s''', method, param_name)
            end
            
            
           
        case 'sigma'
%%
            switch method
                case 'default'
                    if isfield(store, 'selfmerge_theta')
                        theta = store.selfmerge_theta;
                    else
                        theta = D.sumparam(D.scaleparam(w1/(w1+w2),theta1), D.scaleparam(w2/(w1+w2),theta2));
                        store.selfmerge_theta = theta;
                    end
                    value = theta.sigma;

                otherwise
                    error('mvn.selfmerge: Method ''%s'' not recognized for parameter ''%s''', method, param_name)
            end
            
    end
end
