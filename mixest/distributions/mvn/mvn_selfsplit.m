%% |mvn_selfsplit|
% *Note:* This is a private function.
%
% Calculate the initial parameter values for two splitted distributions to
% be substituted for the current distribution. Each parameter should be
% requested by a separate call with an optional specific method.
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

function [value1, value2, store, mixture_store] = ...
    mvn_selfsplit(D, theta, param_name, method, data, store, ...
    mixture_D, mixture_theta, idx, mixture_store) %#ok<INUSL>

    if nargin < 4
        method = 'default';
    end
    if nargin < 6
        store = struct;
    end
    
    switch param_name

        case 'mu'
%%
            switch method
                case 'default'
                    [V,D] = eig(theta.sigma);
                    lambda_max = abs(D(end,end));
                    v_max = V(:,end);
                    mu_add = (sqrt(lambda_max)/2) * v_max;
                    value1 = theta.mu + mu_add;
                    value2 = theta.mu - mu_add;

                otherwise
                    error('mvn.selfsplit: Method ''%s'' not recognized for parameter ''%s''', method, param_name)
            end
            
            
           
        case 'sigma'
%%
            switch method
                case 'default'
                    value1 = theta.sigma / 2;
                    value2 = value1;

                otherwise
                    error('mvn.selfsplit: Method ''%s'' not recognized for parameter ''%s''', method, param_name)
            end
            
    end
end
