function test_plotting
clc
clear
close all

data = 1 + 3.*randn(1, 10000);

D = mvnfactory(1);

options.solver = 'cg';

% options.plotcost = true;
% options.plotgradnorm = true;
options.plotitercount = Inf;
options.plotavg = 5;
figure
options.plotcost.axes = subplot(1,2,1);
options.plotgradnorm.axes = subplot(1,2,2);

options.minibatch.size = 1000;
options.minibatch.iter = 4;

options.crossval = true;
options.crossval.toliter = 10;

options.minibatch.discardhistory = false;

options.miniter = 200;
options.minstepsize = 0;

options.stopfun = @stopfun;
options.verbosity = 2;

options.theta0 = D.randparam();
theta = D.estimate(data, options)

function stopnow = stopfun(D, theta, info, last) %#ok<INUSD>
% pause(0.1)
stopnow = false;
