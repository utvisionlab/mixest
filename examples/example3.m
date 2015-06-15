%% Example 3

%%
% This example visualizes the estimation process of a mixture of von
% Mises-Fisher distributions on sample data.
%

function example3

clc
clear
close all

num = 5;
D = mixturefactory(vmffactory(3), num);

theta = D.randparam();
for k = 1 : num
    theta.D{k}.kappa = theta.D{k}.kappa*10 +10;
end

data = D.sample(theta, 5000);

figure('units', 'normalized', 'outerposition', [0 0 1 1])

% graphical visualization
options.visualization.axes = subplot(2,2,[1 3]);
options.visualization.mixture.colorize = true;
% options.visualization.dataplottype = 'scatter';


% % graphical visualization init
% viz = subplot(2,2,[1 3]);
% plot3(viz, data(1, :), data(2, :), data(3, :), '.', ...
%         'MarkerSize', 4, 'Marker', '.');
% options.stopfun = @(D, theta, info, last) stopfun(D, theta, data, viz);


% plotting options
options.plotcost.axes = subplot(2,2,2);
options.plotgradnorm.axes = subplot(2,2,4);

% main options
options.verbosity = 2;
options.solver = 'default';
% options.crossval = true;
% options.cvtoliter = 100;
options.tolcostdiff = 1e-5;
%options.penalize = true;
%options.minibatchiter = 5;
%options.minibatchdiscardhistory = false;
%options.minibatchsize = 300;
options.regularize = true;
% perform estimation
D.estimate(data, options);


function stopnow = stopfun(D, theta, data, viz)
% update the graphical visualization
colors = hsv(D.num()); % {[ 0 0 1], [0 1 0], [1 0 0]};
weights = D.weighting(theta, data);
[notused, ind] = max(weights,[],1);
hold(viz, 'off');
for k = 1 : D.num()
    mu = theta.D{k}.mu;
    plot3(viz, data(1,ind==k), data(2,ind==k), data(3,ind==k), '.', ...
        'MarkerSize', 4, 'Color', colors(k,:), 'Marker', '.');
    %set(h, 'MarkerSize', 4, 'Color', colors{k}, 'Marker', '.');
    h = line([0 mu(1)], [0 mu(2)], [0 mu(3)], 'Parent', viz, ...
        'LineWidth', 2, 'Color', colors(k,:));
    %set(h, 'LineWidth', 2, 'Color', colors{k}, 'Parent', viz);
    if k == 1
        hold(viz, 'on');
    end
end
drawnow
pause(0.3)
stopnow = false;