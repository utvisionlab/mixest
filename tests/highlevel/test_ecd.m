function test_ecd

n = 10;
nrs = 10000;
x = gammafactory;
D = ecdfactory(n, x);
test_functions(D, nrs)
