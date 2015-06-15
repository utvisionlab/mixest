%% |mvn_visualize|
% *Note:* This is a private function.
%
% von Mises-Fisher distribution visualization
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

function handle_array = vmf_visualize(D, theta, vis_options)

    datadim = D.datadim();
    if ~(datadim==3)
        error('Only 3-D visualization is supported')
    end

    default_options = struct( ...
        'axes', 0, ...
        'color', [0 0 1], ...
        'label', '', ...
        'fixed', false ...
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

    handle_array = [];

    mu = theta.mu;

    h = line([0 mu(1)], [0 mu(2)], [0 mu(3)], 'Parent', ax, ...
        'LineWidth', 2, 'Color', vo.color);
    handle_array = [handle_array; h];

    
    if ~isempty(vo.label)
        h = text(mu(1), mu(2), mu(3), vo.label, 'Parent', ax);
        set(h, ...
            'FontWeight','bold', ...
            'FontSize',24, ...
            'FontUnits','normalized', ...
            'FontName','FixedWidth', ...
...         'BackgroundColor','w', ...
            'HorizontalAlignment','center');
        handle_array = [handle_array; h];
    end

end

