%% |mxe_plot|
% *Note:* This is a private function.
%
% Returns a structure used for plotting costs, etc over iterations
%
% *Syntax*
%
%   mp = mxe_plot(labels)
%   mp = mxe_plot(labels, plot_options)
%   mp = mxe_plot(labels, plot_options, reset_axes)
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

function mp = mxe_plot(labels, plot_options, reset_axes)
% labels is a cell array of strings specifying the legend captions
% (and the number of plots required).
%
% reset_axes is a flag to clear the previous contents of the axes (default
% true). When calling mxe_plot for multiple plots, this Should be set to
% false for all calls except the first one.
%

    if nargin < 3
        reset_axes = true;
    end
    
    ax = plot_options.axes;
    
    % check if the axes has been closed before reaching here (this might
    % occur especially in compound estimations)
    mp.closed = false;
    if ~( ishandle(ax) && strcmp(get(ax,'type'),'axes') )
        mp.closed = true;
    end
    % we need to add the needs_update function
    if ~mp.closed
        
        
        % we accept a single label as a string
        if ischar(labels)
            labels = {labels};
        end

        count = numel(labels);
        mp.h = zeros(count, 1);

        if reset_axes
            cla(ax, 'reset')
        end

        hold(ax, 'all')
        for k = 1 : count
            mp.h(k) = plot(ax, 0, nan, 'DisplayName', labels{k});
        end
        hold(ax, 'off')

        if ischar(plot_options.legend) && strcmp(plot_options.legend, 'default')
            legend(ax, '-DynamicLegend'); % uses the DisplayName properties
        else
            legend(ax, plot_options.legend);
        end

        if ~isempty(plot_options.xlabel)
            xlabel(ax, plot_options.xlabel)
        end
        if ~isempty(plot_options.ylabel)
            ylabel(ax, plot_options.ylabel)
        end
        if ~isempty(plot_options.title)
            title(ax, plot_options.title)
        end
        if plot_options.log
            set(ax, 'Yscale', 'log')
        end
        if plot_options.grid
            grid(ax, 'on')
        end

    end



    mp.needs_update = @needs_update;
    function flag = needs_update(iter, idx)
    % returns false if any of the following conditions are met:
    % * the plot has been closed by user
    % * no update required in current iteration due to options.avgiter
    %
    % idx: index of plot (used when plot count is more than 1)

        if nargin < 2
            idx = 1;
        end
        
        flag = true;
        
        if mp.closed || ~ishandle(mp.h(idx))
            flag = false;
            return
        end
        
        if iter==0 || mod(iter, plot_options.avgiter)~=0
            flag = false;
            return
        end
    end

    mp.update = @update;
    function update(xData, yData, idx)
    % updates a plot with xData, yData
    % idx: index of plot (used when plot count is more than 1)
    
        if nargin < 3
            idx = 1;
        end
        set(mp.h(idx), 'XData', xData, 'YData', yData);
    end

end