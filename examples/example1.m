%% Example 1

%%
% This example visualizes the estimation process of a mixture of three
% Gaussians on sample 2-D data.
% 

function example1

clc
clear
close all

load data2d

num = 3;
D = mixturefactory(mvn2factory(2), num);


% graphical visualization
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1])
options.visualization.axes = subplot(2,2,[1 3]);
options.visualization.mixture.colorize = true;
% options.visualization.dataPlotType = 'scatter';
options.visualization.stopOnClose = true;

% plotting
options.plotCost.axes = subplot(2,2,2);
options.plotGradNorm.axes = subplot(2,2,4);

% main options
options.verbosity = 1;
options.solver = 'lbfgs';
options.sgd.stepsize = 4;
% options.crossVal = true;
% options.crossVal.tolIter = 100;
% options.tolCostDiff = -Inf;
% options.penalize = true;
% options.minibatch.iter = 5;
% options.minibatch.discardHistory = false;
% options.minibatch.size = 300;
% options.regularize = true;

% perform estimation
D.estimate(data, options)
