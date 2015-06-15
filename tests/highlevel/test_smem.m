function test_smem

clc
clear
close all

f = load('data2d');
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

% % SMEM options
% options.sm.maxcands = 5;
% options.sm.ComponentD = mvnfactory2(2);

% % inner estimation options
% options.inner.full.maxIter = 500;
% options.inner.full.tolCostDiff = 1e-7;
% options.inner.partial.maxIter = 500;
% options.inner.partial.tolCostDiff = 1e-7;


% run
smem(data, 3, options)
