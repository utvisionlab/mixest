function test_crossvalidation

data = 1 + 3.*randn(1, 1000);

D = mvnfactory(1);

options.solver = 'cg';
options.crossval = true;

theta = D.estimate(data, options)

assertElementsAlmostEqual(theta.mu, 1, 'absolute', 0.5)
assertElementsAlmostEqual(theta.sigma, 9, 'absolute', 0.5)
