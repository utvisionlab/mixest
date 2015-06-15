function test_factorial

n = 10;
nrs = 1000;
k = 2;
s = mixturefactory(mvnfactory(1), k);
D = factorialfactory(s, n);

test_functions(D, nrs);