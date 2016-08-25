%% |mixture_mergeinit|
% *Note:* This is a private function.
%
% Calculate the initial parameters for a merged mixture component to be
% substituted for two given components
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

function [theta, store] = mixture_mergeinit(D, idx1, idx2, theta, options, data, store)

    if nargin < 7
        store = struct;
    end
    if nargin>=5 && isfield(options, 'sm')
        if isa(options.sm.mergeinit, 'function_handle')
            [theta, store] = options.sm.mergeinit(D, idx1, idx2, theta, options, data, store);
            return
        end
        options_mergeinit = options.sm.mergeinit;
    else        
        options_mergeinit = 'default';
    end
    
%%

    % component parameters
    Didx1 = D.component(idx1);
    Didx1Name = Didx1.name();
    param_names = fieldnames(theta.D{idx1});
    
    if isfield(store, 'componentStores')
        cmptStore = store.componentStores{idx};
    else
        cmptStore = struct;
    end
    
    if isfield(Didx1, 'selfmerge')
        for m = 1 : numel(param_names)
            param_name = param_names{m};
            method = get_method(options_mergeinit, Didx1Name, param_name);

            [value, cmptStore, store] = ...
                Didx1.selfmerge(theta.D{idx1}, theta.D{idx2}, param_name, ...
                theta.p(idx1), theta.p(idx2), method, data, cmptStore, ...
                D, theta, idx1, idx2, store);
            theta.D{idx1}.(param_name) = value;
        end
        
    else
        % if the component distribution does not define a 'selfmerge'
        % function, we use the default merging method
        method = get_method(options_mergeinit, Didx1Name, '');
        if strcmp(method, 'default')
            if true
                options.verbosity = 0;
                DidxM = D.component(idx1);
                % penalizer_theta
                if options.penalize
                    if isempty(options.penalizertheta)
                        options.penalizertheta = DidxM.penalizerparam(data);
                    end
                end
                data = mxe_readdata(data);
                h = D.weighting(theta, data);
                data.weight = sum(h([idx1 idx2],:),1);
                options.theta0 = theta.D{idx1};
                theta.D{idx1} = DidxM.estimate(data, options);
            else
                w1 = theta.p(idx1);
                w2 = theta.p(idx2);
                theta.D{idx1} = Didx1.sumparam(...
                    Didx1.scaleparam(w1/(w1+w2),theta.D{idx1}), ...
                    Didx1.scaleparam(w2/(w1+w2),theta.D{idx2}));
            end
        else
            error('mixture.mergeinit: Method ''%s'' not recognized for component(%d) parameters', method, idx)
        end
    end
        
    theta.D(idx2) = [];
    
%%

    % component weights (Note: this should be done after component parameters since we use theta.p there)
    param_name = 'p';
    method = get_method(options_mergeinit, D.name(), param_name);
    switch method
        % every method should set p for the weight of the merged component
        case 'default'
            p = theta.p(idx1) + theta.p(idx2);
            
        otherwise
            error('mixture.mergeinit: Method ''%s'' not recognized for parameter ''%s''', method, param_name)
    end
    theta.p(idx1) = p;
    theta.p(idx2) = [];
    
end



function method = get_method(options_mergeinit, DName, param_name)
% extracts the method from options_mergeinit

    method = options_mergeinit;
    if isfield(method, DName) % isfield checks for isstruct also
        method = method.(DName);
        if isfield(method, param_name)
            method = method.(param_name);
        else
            method = 'default';
        end
    else
        method = 'default';
    end
end
