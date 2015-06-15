%% |mvn_visualize|
% *Note:* This is a private function.
%
% Multi-variate-normal distribution visualization
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

function handle_array = mvn_visualize(D, theta, vis_options)

    datadim = D.datadim();
    if ~(datadim==2 || datadim==3)
        error('Only 2-D and 3-D visualization is supported')
    end

    default_options = struct( ...
        'axes', 0, ...
        'color', [0 0 1], ...
        'label', '', ...
        'showcenter', true, ...
        'ellipse_n', 100, ...
        'ellipsoid_n', 20, ...
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

    sigma = theta.sigma;
    mu = theta.mu;

    [eigV, eigD] = eig(sigma);
    lambda = sqrt(diag(eigD));
    
    % draw cluster ellipse
    if datadim >= 3
        h = drawEllipsoid(ax, mu, eigV, lambda, vo.ellipsoid_n);
        set(h, 'EdgeColor', vo.color, 'LineWidth', 2);
    else
        h = drawEllipse(ax, mu, eigV, lambda, vo.ellipse_n);
        set(h, 'Color', vo.color, 'LineWidth', 2);
    end
    if vo.fixed
        set(h, 'LineStyle', '--');
    end
    handle_array = [handle_array; h];

    
    
    % draw eigenvectors
    if vo.showcenter
        for i = 1:datadim
            Vi = eigV(:,i);
            Di = lambda(i);

            X1 = mu - Di.*Vi;
            X2 = mu + Di.*Vi;
            if datadim >= 3
                h = line([X1(1) X2(1)], [X1(2) X2(2)], [X1(3) X2(3)], 'Parent', ax);
            else
                h = line([X1(1) X2(1)], [X1(2) X2(2)], 'Parent', ax);
            end
            set(h, 'Color', vo.color, 'LineWidth', 2);
            if vo.fixed
                set(h, 'LineStyle', '--');
            end
            handle_array = [handle_array; h]; %#ok<AGROW>
        end
    end

    if ~isempty(vo.label)
        if datadim >= 3
            h = text(mu(1), mu(2), mu(3), vo.label, 'Parent', ax);
        else
            h = text(mu(1), mu(2), vo.label, 'Parent', ax);
        end
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



function h = drawEllipse(ax, mu, eigV, lambda, ellipse_n)

    V2 = eigV(:,2);
    angle = atan2(V2(2), V2(1));
    if(angle < 0)
        angle = angle + 2*pi;
    end
    a = 3 * lambda(2);
    b = 3 * lambda(1);

    X = calcEllipse(mu, a, b, angle, ellipse_n);

    h = plot(ax, X(:,1), X(:,2));
end
function X = calcEllipse(mu, a, b, angle, ellipse_n)

    % find coordinates on the canonical ellipse
    phi_vec = transpose(linspace(0, 2*pi, ellipse_n));
    x  = a * cos(phi_vec);
    y  = b * sin(phi_vec);

    % apply rotation
    R = [cos(angle) sin(angle); -sin(angle) cos(angle)];
    X = [x y] * R;
    
    % apply translation
    X = X + repmat(mu(:).', [size(X,1),1]);
end



function h = drawEllipsoid(ax, mu, eigV, lambda, ellipsoid_n)

    [x, y, z] = calcEllipsoid(mu, eigV, 3 * lambda, ellipsoid_n);
    
    h = surf(ax, x, y, z, 'FaceColor', 'none');
end
function [x, y, z] = calcEllipsoid(mu, eigV, lambda, ellipsoid_n)

    % find coordinates on the canonical ellipsoid
    [x, y, z] = ellipsoid(0,0,0,lambda(1),lambda(2),lambda(3),ellipsoid_n);

    % apply rotation and translation
    X = kron(eigV(:,1), x) + kron(eigV(:,2), y) + kron(eigV(:,3), z);
    n = size(X, 2);
    x = X(1:n, :) + mu(1); 
    y = X(n+1:2*n, :) + mu(2); 
    z = X(2*n+1:end, :) + mu(3);
end


