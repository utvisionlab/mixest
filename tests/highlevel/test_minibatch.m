function test_minibatch

data = 1 + 3.*randn(1, 1000);

D = mvnfactory(1);

options.solver = 'cg';
options.minibatch.size = 10;
options.minibatch.discardhistory = true;
options.minibatch.iter = 1;

options.crossval = true;

theta = D.estimate(data, options)
