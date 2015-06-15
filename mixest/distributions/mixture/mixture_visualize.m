%% |mixture_visualize|
% *Note:* This is a private function.
%
% Mixture distribution visualization
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

function handle_array = mixture_visualize(D, theta, vis_options)

    default_options = struct( ...
        'axes', 0, ...
        'colorize', false, ...
        'colormap', [], ...
        'showlabels', true, ...
        'data', [], ...
        'hdata', 0, ...
        'dataplottype', 'patch' ... % 'patch' (faster), 'scatter' (slower)
        );
    
    if nargin < 3
        vo = default_options;
    else
        vo = mxe_setfields(default_options, vis_options);
    end
    if ~( ishandle(vo.axes) && strcmp(get(vo.axes,'type'),'axes') )
        figure
        vo.axes = gca;
    end
    ax = vo.axes;
    
    if ~isempty(vo.data)
        if vo.hdata == 0
            data = mxe_readdata(vo.data);
            datadim = data.dim;
            data = data.data;
            %TODO redundant
            switch vo.dataplottype
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
            hold(ax, 'on')
            vo.hdata = h;
            clear data
        end
    end
    
    handle_array = [];

    num = D.num();
    numtotal = D.numtotal();
    if numtotal > num
        numWeighting = num + 1; %TODO weighting returns this number of rows
    else
        numWeighting = num;
    end

    % update colors
    if vo.colorize
        cmap = vo.colormap;
        if isempty(cmap)
            cmap = hsv(numWeighting); % cluster color map (each row contains the color related to a cluster)
        elseif size(cmap,1) < numWeighting
            cmap = repmat(cmap, numWeighting, 1);
            cmap = cmap(1:numWeighting,:);
        end
        
        if ~isempty(vo.data)
            hX = D.weighting(theta, vo.data);
            R = hX.' * cmap(:, 1);
            G = hX.' * cmap(:, 2);
            B = hX.' * cmap(:, 3);
            colors = [R G B];
            
            switch vo.dataplottype
                case 'patch'
                    set(vo.hdata, 'FaceVertexCData',colors, 'FaceColor','flat')
                case 'scatter'
                    set(vo.hdata, 'CData', colors)
            end
        end
    end

    
    % cluster visualization
    for j = 1:numtotal
        
        if j > num
            Dj = D.fixedD(j-num);
            fixedtheta = D.fixedparam();
            theta_j = fixedtheta.D{j-num};
            fixed = true;
        else
            Dj = D.varD(j);
            theta_j = theta.D{j};
            fixed = false;
        end
        
        DjName = Dj.name(); % custom options for components may be passed through vo.(DjName)
        
        vo.(DjName).axes = ax;
        vo.(DjName).fixed = fixed;
        
        if vo.colorize
            if j > num
                vo.(DjName).color = cmap(num+1,:);
            else
                vo.(DjName).color = cmap(j,:);
            end
        end
        
        if vo.showlabels
            vo.(DjName).label = num2str(j);
        end
        
        if isfield(Dj, 'visualize')
            h = Dj.visualize(theta_j, vo.(DjName));
            handle_array = [handle_array; h]; %#ok<AGROW>
        else
            warning('visualize not implemented for distribution %s (component %d)', DjName, j)
        end
    end

end
