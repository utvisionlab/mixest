%% |mxe_inneroptions|
% *Note:* This is a private function.
%
% Read options for inner estimations
%
% *Syntax*
%
%   inner_options = mxe_inneroptions(options, inner_defaults, subfield)
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

function inner_options = mxe_inneroptions(options, inner_defaults, subfield)
% |options|: full options, possibly containing options.inner fields given by the user
% |inner_defaults| (optional): specific defaults for the inner options
% |subfield| (optional): specific sub-field of options.inner, used to address specific inner estimation functions between various functions
%
% The estimation functions that call other estimation functions inside,
% should use this function to get the options for the inner estimation
% functions. They can give specific defaults using the |inner_defaults|
% argument and also they may give |subfield|, defining a specific name
% distinguishing the inner estimation from any other inner estimations, in
% order to enable the user to change the options for that specific inner
% estimation by |options.inner.(subfield).(option)|.
%
% The priorities for setting the output options are as follows:
% 1. |options.inner| (if set by the user)
% 2. |inner_defaults|
% 3. |options| (outer options)

    
    if isfield(options, 'inner')
        user_options = options.inner;
        options = rmfield(options, 'inner');
        
        if nargin>2 && isfield(user_options, subfield)
            user_options = user_options.(subfield);
        end
    else
        user_options = [];
    end
    
    if nargin > 1
        options = mxe_setfields(options, inner_defaults);
    end
    
    inner_options = mxe_setfields(options, user_options);

end