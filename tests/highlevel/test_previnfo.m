function test_previnfo

data = randn(1,1000)*2+3;
D = mvnfactory(1);

% get some nice output
options.solver = 'cg';
options.plotcost.axes = gca;
options.plotitercount = inf;
options.verbosity = 2;
options.maxiter = 50;
options.minibatch.size = 10;

% don't stop too soon!
options.tolgradnorm = -inf;
options.tolcost = -inf;
options.tolcostdiff = -inf;
options.minstepsize = -inf;

% crossval
options.crossval = true;
options.miniter = 100;

% gradnorm
options.plotgradnorm.axes = options.plotcost.axes;

% first estimation
[theta, D, info, options] = D.estimate(data, options);

% test previnfo on second estimation
options.previnfo = info;
options.theta0 = theta;
options.maxiter = 100;
[theta, D, info] = D.estimate(data, options);
