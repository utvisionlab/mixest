%% |mxe_setfields|
% *Note:* This is a private function.
%
% Sets options structure recursively, converting the names to lower case
%
% *Syntax*
%
%   options = mxe_setfields(options, given_options)
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

function options = mxe_setfields(options, given_options)
% sets the given fields from given_options in options.
% struct-type options are handled recursively.
% given field names are converted to lower case.

    if isempty(given_options)
        return
    end
    
    % MATLAB R14: Assigning Nonstructure Variables As Structures Displays Warning
    if ~isstruct(options)
        options = struct;
    end
    
    given_fields = fieldnames(given_options);
    for n = 1:numel(given_fields)
        fn = given_fields{n};
        lfn = lower(fn);
        
        issuboption = isstruct(given_options.(fn));
        if issuboption && numel(given_options.(fn)) > 1 % for info structs passed as options (previnfo)
            issuboption = false;
        end
        if issuboption && isfield(given_options.(fn), 'name') % for distribution structs passed as options (we shouldn't make their field names lower case)
            issuboption = false;
        end
        
        if issuboption
            if ~isfield(options, lfn)
                options.(lfn) = struct;
            end
            options.(lfn) = mxe_setfields(options.(lfn), given_options.(fn));
        else
            options.(lfn) = given_options.(fn);
        end
    end
end

