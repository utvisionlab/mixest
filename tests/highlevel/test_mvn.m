function test_mvn

n = 10;
nrs = 10000;
D = mvnfactory(n);

test_functions(D, nrs)

% Calculating kl-divergence
theta = D.randparam();
theta2 = D.randparam();
D.kl(theta, theta2);