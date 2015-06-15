function test_cem_spiral

clc
clear
close all

f = load('spiral');
data = f.data;

% visualization and plotting options
figure('Units', 'normalized', 'OuterPosition', [0 0.3 1 0.7])
options.visualization.axes = subplot(2,2,[1 3]);
options.plotCost.axes = subplot(2,2,2);
options.plotGradNorm.axes = subplot(2,2,4);
options.visualization.mixture.colorize = true;
options.visualization.stopOnClose = true;

% common estimation options
options.verbosity = 1;
options.solver = 'lbfgs';
options.tolCostDiff = 1e-3;

% split-and-merge options
% options.sm.splitCriterion = 'kl';
% options.sm.mergeCriterion = 'kl';
% options.sm.splitInit = 'default';
% options.sm.mergeInit = 'default';
% options.sm.numInit = 20;
options.sm.numMin = 1;
options.sm.numMax = 40;
options.sm.mergeMaxCands = 12;
options.sm.splitMaxCands = 12;
options.sm.maxFail = 2;
% options.sm.ComponentD = mvnfactory2(2);

% inner estimation options
options.inner.full.maxIter = 500;
options.inner.full.tolCostDiff = 1e-7;
options.inner.partial.maxIter = 500;
options.inner.partial.tolCostDiff = 1e-7;


% run
cem(data, options)
