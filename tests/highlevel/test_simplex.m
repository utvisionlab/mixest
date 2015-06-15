function test_simplex
profile on
q= 10;
M = simplexKfactory(q,2);

problem.M = M;
disp('Checking ll gradient using egrad2rgrad');
problem.costgrad = @(x)costgrad(M,x);
figure(1),checkgradient(problem)

profile viewer
profile off
function [lik,dll] = costgrad(M,x)
lik = sum(x.^3);
dll = 3*x.^2;
dll = M.egrad2rgrad(x,dll);