function test_mixture

n = 10;
nrs = 10000;
x = mvnfactory(n);
D = mixturefactory(x, 3);
test_functions(D, nrs)

% Calculating kl-divergence
D.kl(theta);
