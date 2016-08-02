%% |mixture_splitinit|
% *Note:* This is a private function.
%
% Calculate the initial parameters for two splitted mixture components to
% be substituted for the given component
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

function [theta, store] = mixture_splitinit(D, idx, theta, options, data, store)

    if nargin < 6
        store = struct;
    end
    if nargin>=4 && isfield(options, 'sm')
        if isa(options.sm.splitinit, 'function_handle')
            [theta, store] = options.sm.splitinit(D, idx, theta, options, data, store);
            return
        end
        options_splitinit = options.sm.splitinit;
    else        
        options_splitinit = 'default';
    end
    num = D.num();
    numfixed = D.numfixed();
    
%%

    % component parameters
    Didx = D.component(idx);
    DidxName = Didx.name();
    param_names = fieldnames(theta.D{idx});
    theta.D = [theta.D(:); theta.D{idx}];
    
    if isfield(store, 'componentStores')
        cmptStore = store.componentStores{idx};
    else
        cmptStore = struct;
    end
    
    if isfield(Didx, 'selfsplit')
        for m = 1 : numel(param_names)
            param_name = param_names{m};
            method = get_method(options_splitinit, DidxName, param_name);

            [value1, value2, cmptStore, store] = ...
                Didx.selfsplit(theta.D{idx}, param_name, method, data, cmptStore, D, theta, idx, store);

            theta.D{idx}.(param_name) = value1;
            theta.D{end}.(param_name) = value2;
        end
        
    else
        % if the component distribution does not define a 'selfsplit'
        % function, we add random values to each parameter
        if true
            DidxM = Didx.M;
            DidxRandvec = DidxM.randvec(theta.D{idx});
            theta.D{end} = DidxM.retr(theta.D{idx},DidxRandvec, 2e-10);
            theta.D{idx} = DidxM.retr(theta.D{idx},DidxRandvec, 1e-10);
        else
            for m = 1 : numel(param_names)
                param_name = param_names{m};
                method = get_method(options_splitinit, DidxName, param_name);
                
                if strcmp(method, 'default')
                    oldvalue = theta.D{idx}.(param_name);
                    value1 = oldvalue + rand(size(oldvalue)) * norm(oldvalue); %TODO what amount of noise is best?
                    value2 = oldvalue + rand(size(oldvalue)) * norm(oldvalue);
                    
                    theta.D{idx}.(param_name) = value1;
                    theta.D{end}.(param_name) = value2;
                else
                    error('mixture.splitinit: Method ''%s'' not recognized for component(%d) parameter ''%s''', method, idx, param_name)
                end
            end
        end
    end
    
%%

    % component weights (Note: this should be done after component parameters since we use theta.p there)
    param_name = 'p';
    method = get_method(options_splitinit, D.name(), param_name);
    switch method
        % every method should set p1,p2 for the weights of the first and
        % second splitted components respectively
        case 'default'
            p1 = theta.p(idx) / 2;
            p2 = p1;
            
        otherwise
            error('mixture.splitinit: Method ''%s'' not recognized for parameter ''%s''', method, param_name)
    end
    if numfixed > 0
        w = theta.p(num+1);
        theta.p(idx) = p1;
        theta.p = [theta.p(1:num); p2; w];
    else
        theta.p(idx) = p1;
        theta.p = [theta.p(:); p2];
    end

end



function method = get_method(options_splitinit, DName, param_name)
% extracts the method from options_splitinit

    method = options_splitinit;
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
