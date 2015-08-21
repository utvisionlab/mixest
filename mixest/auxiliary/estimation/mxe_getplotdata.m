%% |mxe_getplotdata|
% *Note:* This is a private function.
%
% Extract plot data from info
%
% *Syntax*
%
%   [ploty, plotx] = mxe_getplotdata(fieldname, plot_options, info, last)
%   [ploty, plotx] = mxe_getplotdata(fieldname, plot_options, info, last, previter, allinfo, allinfolast)
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

function [ploty, plotx] = mxe_getplotdata(fieldname, plot_options, info, last, previter, allinfo, allinfolast)

    if nargin < 5
        

        if plot_options.itercount <= last
            % available info is enough
            plotx = [info(last-plot_options.itercount+1 : last).iter];
            ploty = [info(last-plot_options.itercount+1 : last).(fieldname)];
        else
            % available info is less than iterCount
            plotx = [info(1:last).iter];
            ploty = [info(1:last).(fieldname)];
        end
        
        
    else
        
        
        iterCount = plot_options.itercount;
        if iterCount == Inf
            % iterCount==Inf => plot all info
            if previter > 0
                plotx = [allinfo.iter, info(2:last).iter];
                ploty = [allinfo.(fieldname), info(2:last).(fieldname)];
            else
                plotx = [info(1:last).iter];
                ploty = [info(1:last).(fieldname)];
            end
        elseif iterCount<last || ...
                (previter==0 && iterCount==last)
            % all info to be plotted is in current info structure array
            plotx = [info(last-iterCount+1 : last).iter];
            ploty = [info(last-iterCount+1 : last).(fieldname)];
        elseif previter > 0
            % we need info from previous runs
            startidx = max(allinfolast - (iterCount - last), 1);
            plotx = [allinfo(startidx:allinfolast).iter, info(2:last).iter];
            ploty = [allinfo(startidx:allinfolast).(fieldname), info(2:last).(fieldname)];
        else
            % available info is less than iterCount
            plotx = [info(1:last).iter];
            ploty = [info(1:last).(fieldname)];
        end
        
        
    end
    

    % averaged plotting
    if plot_options.avgiter > 1
        [plotx, ploty] = mxe_plotavg(plotx, ploty, plot_options.avgiter);
    end
end
