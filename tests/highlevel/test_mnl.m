function test_mnl


D = mnlfactory(2,2);
data = rand(2,100);
cat = ceil(rand(1,100)*2);
test_functions(D, [data;cat])

% Calculating kl-divergence
theta = D.randparam();
theta2 = D.randparam();
D.kl(theta, theta2);