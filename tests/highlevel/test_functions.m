function thetaO = test_functions(D, nrs)

profile on

% Generating random samples
theta = D.randparam();

if isscalar(nrs)
    data = D.sample(theta, nrs);
else
    data = nrs;
end
%figure(5), plot(data(1,:),data(2,:),'.');

if true
    % Checking penalizer
    disp('Checking penalizer');
    Pparam = D.penalizerparam(data);
    problem.M = D.M;
    problem.costgrad = @(x)costgradpen(D,x,Pparam);
    figure(1),checkgradient(problem)
end

if false
    disp('Checking ll gradient without egrad2rgrad');
    problem.M = D.M;
    problem.costgrad = @(x) D.ll(x, data);
    figure(1), checkgradient(problem)
end

if false
    % Checking gradient with respect to data
    disp('Checking ll gradient w.r.t data');
    problem.M = euclideanfactory(size(data,1));
    problem.costgrad = @(data)costgraddata(D,theta,data);
    figure(2),checkgradient(problem)
end

if false
    % Checking llgrad
    disp('Checking ll gradient using egrad2rgrad');
    problem.M = D.M;
    problem.costgrad = @(x)costgrad(D,x,data);
    figure(3),checkgradient(problem)
end

if true
    % Checking llgrad
    disp('Checking ll plus pen gradient using egrad2rgrad');
    problem.M = D.M;
    problem.costgrad = @(x)costgradupen(D,x,data,Pparam);
    figure(5),checkgradient(problem)
end
% Checking weighted llgrad
disp('Checking weighted ll gradient using egrad2rgrad');
problem.M = D.M;
problem.costgrad = @(x)costgrad(D,x,struct('data',data,'weight',rand(1,nrs)));
figure(4),checkgradient(problem)

% setting the optimization procedure to lbfgs
options.verbosity = 2;
options.solver = 'lbfgs';
options.tolcost = -Inf;
options.tolcostdiff = 1e-4;
options.penalize = true;
thetaO = D.estimate(data, options);

if false
    % Calculating entropy
    disp(['Theoretical Entropy = ' num2str(D.entropy(theta))]);
    disp(['Numerical Entropy = ' num2str(-D.ll(theta,data)/nrs)]);
    profile viewer
    profile off
end

function [lik,dld] = costgraddata(d,x,data)
[lik,store] = d.ll(x,data);
dld = d.llgraddata(x,data,store);

function [lik,dll] = costgradpen(d,x,ppar)
[lik,store] = d.penalizercost(x,ppar);
dll = d.penalizergrad(x,ppar,store);
dll = d.M.egrad2rgrad(x,dll);

function [lik,dll] = costgrad(d,x,data)
[lik,store] = d.ll(x,data);
dll = d.llgrad(x,data,store);
dll = d.M.egrad2rgrad(x,dll);


function [lik,dll] = costgradupen(d,x,data,ppar)
options = mxe_options;
options.penalize = true;
options.penalizertheta = ppar;
[lik,dll] = mxe_costgrad(d,x,data,options);