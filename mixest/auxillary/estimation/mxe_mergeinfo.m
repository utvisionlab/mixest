%% |mxe_mergeinfo|
% *Note:* This is a private function.
%
% Merge info structure arrays, considering memory pre-allocation
%
% *Syntax*
%
%   [info, last] = mxe_mergeinfo(info, last, newinfo)
%   [info, last] = mxe_mergeinfo(info, last, newinfo, maxiter)
%   [info, last] = mxe_mergeinfo(info, last, newinfo, maxiter, field_names_differ)
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

function [info, last] = mxe_mergeinfo(info, last, newinfo, maxiter, field_names_differ)
% maxiter: max size for info
% field_names_differ: flag to use the hard way

    % if info is empty, return the newinfo
    if isempty(info)
        info = newinfo;
        last = numel(newinfo);
        return
    end
    
    % if newinfo is empty, return the info untouched
    if isempty(newinfo)
        return
    end
    
    if nargin < 5
        field_names_differ = false;
    end
    if nargin < 4
        maxiter = Inf;
    end
    
    % update allinfo and pre-allocate if necessary
    allinfosize = numel(info);
    infosize = numel(newinfo)-1;
    tobelast = last + infosize;
    if tobelast > allinfosize
        allinfosize = min(tobelast+10000, maxiter+1);
        info(allinfosize).iter = [];
    end
    
    if ~field_names_differ
        % perform the merge
        info(last+1 : tobelast) = newinfo(2:end);
        last = tobelast;
        return
    end
    
    % manage different field names
    fn_info = fieldnames(info);
    fn_newinfo = fieldnames(newinfo);
    
    % find field names in newinfo that are not in info
    idx_newinfo = find(~ismember(fn_newinfo, fn_info));
    % add them to info
    for m = idx_newinfo
        info(1).(fn_newinfo{m}) = [];
    end
    
    % find field names in info that are not in newinfo
    idx_info = find(~ismember(fn_info, fn_newinfo));
    % add them to newinfo
    for m = idx_info
        newinfo(1).(fn_info{m}) = [];
    end
    
    % we need to make the struct fields in the same order
    info = orderfields(info);
    newinfo = orderfields(newinfo);

    % perform the merge
    info(last+1 : tobelast) = newinfo(2:end);
    last = tobelast;
end
