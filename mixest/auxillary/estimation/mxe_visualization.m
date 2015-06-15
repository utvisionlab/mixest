%% |mxe_visualization|
% *Note:* This is a private function.
%
% Returns a structure used for visualizing the estimation process
%
% *Syntax*
%
%   vis = mxe_visualization(data, vis_options)
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

function vis = mxe_visualization(data, vis_options)

    data = mxe_readdata(data);
    datadim = data.dim;
    data = data.data;

    ax = vis_options.axes;
    
    % check if the axes has been closed before reaching here (this might
    % occur especially in compound estimations)
    % Note: we need to add the update function
    if ~closed()
        
        %TODO don't re-draw the data in compound estimations
        cla(ax, 'reset')

        % visualize data
        switch vis_options.dataplottype
            case 'patch'
                if datadim >= 3
                    Edge = .02; %TODO
                    X=[data(1,:); data(1,:)+Edge; data(1,:)+Edge; data(1,:)];
                    Y=[data(2,:); data(2,:); data(2,:)+Edge; data(2,:)+Edge];
                    Z=[data(3,:); data(3,:); data(3,:); data(3,:)];
                    h = patch(X,Y,Z,'g', 'Parent',ax, 'EdgeColor','none');  
                    view(ax, 3)
                else
                    edge = .1;
                    X=[data(1,:); data(1,:)+edge; data(1,:)+edge; data(1,:)];
                    Y=[data(2,:); data(2,:); data(2,:)+edge; data(2,:)+edge];
                    h = patch(X,Y,'g', 'Parent',ax, 'EdgeColor','none');  
                end
                
            case 'scatter'
                if datadim >= 3
                    h = scatter3(ax, data(1,:), data(2,:), data(3,:), 30, 'g', 'fill');
                else
                    h = scatter(ax, data(1,:), data(2,:), 30, 'g', 'fill');
                end
        end
        drawnow
        hold(ax, 'on')
        handle_array = [];
    end
    

    vis.update = @update;
    function stop = update(theta, D)

        stop = false;
        
        if vis.closed() || ~ishandle(h)
            if vis_options.stoponclose
                stop = true;
            end
            return
        end

        Dname = D.name();
        if isfield(vis_options, Dname)
            vo = mxe_setfields(vis_options, vis_options.(Dname));
        else
            vo = vis_options;
        end
        
        vo.data = data;
        vo.hdata = h;
        
        axes(ax)
        delete(handle_array)
        hold(ax, 'on')
        handle_array = D.visualize(theta, vo);
        hold(ax, 'off')

        drawnow
    end

    vis.closed = @closed;
    function flag = closed()
        flag = ~( ishandle(ax) && strcmp(get(ax,'type'),'axes') );
    end
    
end
