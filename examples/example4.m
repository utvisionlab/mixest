%% Example 4

%%
% This example visualizes the estimation process of a modified version of
% the Competitive EM (CEM) algorithm (Zhang et al., 2004) on sample 2-D
% data.
% 

function example4

clc
clear
close all

f = load('data_sm');
data = f.data;

% visualization and plotting options
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1])
options.visualization.axes = subplot(2,2,[1 3]);
options.plotCost.axes = subplot(2,2,2);
options.plotGradNorm.axes = subplot(2,2,4);
options.visualization.mixture.colorize = true;
options.visualization.stopOnClose = true;

% common estimation options
options.verbosity = 1;
options.solver = 'lbfgs';
options.tolCostDiff = 1e-3;

% % split-and-merge options
% options.sm.splitCriterion = 'kl';
% options.sm.mergeCriterion = 'kl';
% options.sm.splitInit = 'default';
% options.sm.mergeInit = 'default';
% options.sm.numInit = 1;
% options.sm.numMin = 1;
% options.sm.numMax = 15;
% options.sm.mergeMaxCands = 5;
% options.sm.splitMaxCands = 5;
% options.sm.maxFail = 2;
% options.sm.ComponentD = mvnfactory2(2);


% run
cem(data, options)
