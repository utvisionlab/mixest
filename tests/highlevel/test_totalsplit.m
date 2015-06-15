function theta = test_totalsplit

clc
clear
close all
warning off all
%reset(RandStream.getGlobalStream);

load data2d

figure('units', 'normalized', 'outerposition', [0 0.3 1 0.7])
% graphical visualization
options.visualization.axes = subplot(2,2,[1 3]);
% options.visualization.mixture.colorize = true;
options.visualization.stopOnClose = true;
% plotting options
options.plotCost.axes = subplot(2,2,2);
options.plotGradNorm.axes = subplot(2,2,4);


% main options
options.verbosity = 1;
options.debug = true;

% options.solver = 'cg';
% options.crossVal = true;
% options.crossVal.tolIter = 100;
options.regularize = true;
options.penalize = true;
options.miniter = 20;
%options.minibatch.iter = 15;
options.tolCostDiff = 1e-4;%size(data,2) * 1e-4;
% options.tolCostDiff = -Inf;
% options.minibatch.discardHistory = false;
%options.minibatch.size = 500;
 options.sm.costType = 2;
% options.inner.full.maxIter = 10;
% options.inner.partial.maxIter = 2;

options.sm.splitcriterion ='mll';
options.sm.mergecriterion = 'overlap';

options.sm.ComponentD = mvn2factory(2);

target_num = 4;
theta = totalsplit(data, target_num, options);

totalmerge(data, theta, 1, options);
