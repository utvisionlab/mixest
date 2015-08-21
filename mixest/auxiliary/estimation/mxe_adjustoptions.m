%% |mxe_adjustoptions|
% *Note:* This is a private function.
%
% Makes the required adjustments to |options| fields like |maxiter| when
% using |options.previnfo| to continue from a previous estimation.
%
% *Syntax*
%
%   options = mxe_adjustoptions(options, stats)
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

function options = mxe_adjustoptions(options, stats)
% stats is normally options.previnfo(end)

    options.miniter = options.miniter + stats.iter;
    options.maxiter = options.maxiter + stats.iter;
    if isfield(stats, 'time')
        options.maxtime = options.maxtime + stats.time;
    end
    if isfield(stats, 'costevals')
        options.maxcostevals = options.maxcostevals + stats.costevals;
    end

end