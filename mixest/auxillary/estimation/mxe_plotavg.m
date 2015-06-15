%% |mxe_plotavg|
% *Note:* This is a private function.
%
% Function used in averaged plotting
%
% *Syntax*
%
%   [plotx, ploty] = mxe_plotavg(plotx, ploty, avgiter)
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

function [plotx, ploty] = mxe_plotavg(plotx, ploty, avgiter)

    % make the number of elements of plot data a multiple of
    % plotavg, by removing extra elements from its start.
    n = numel(plotx);
    m = mod(n, avgiter);
    if m > 0
        n = n - m;
        plotx = plotx(:, end-n+1:end);
        ploty = ploty(:, end-n+1:end);
    end

    % average every avgiter iterations
    temp = reshape(ploty, avgiter, []);
    ploty = mean(temp, 1);

    % update plotx also
    idx = avgiter : avgiter : n;
    plotx = plotx(:, idx);
end