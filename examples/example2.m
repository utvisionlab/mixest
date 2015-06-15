%% Example 2

%%
% This example uses a mixture of multinomial logistic experts to classify
% Fisher's iris data.
%

function example2

clc
clear
%close all

load data_iris

% random permute the iris data
index = randperm(150);
% load rndprm150
data = data(:, index);

% Use 120 data for training
data_train = data(:, 1:120);
data_test = data(:, 121:150);
datadim = size(data_train,1) - 1;
num = max(data(end,:),[], 2);

% create a mixture of multivariate logit
nummix = 3; % number of components in gate
E = mnlfactory(datadim, num); % expert
G = softmaxfactory(nummix, datadim); % gate
D = moefactory(G, E);
%D = E;
% figure('units', 'normalized', 'outerposition', [0 0 1 1])

% plotting options
options.plotcost = true;


% main options
options.verbosity = 2;
options.solver = 'lbfgs';
options.tolcostdiff = 1e-6;
options.crossval = true;
options.minibatch.size = 10;
options.maxiter = 200;
options.minibatch.discardhistory = false;
options.crossval.toliter = 100;

% perform estimation
theta = D.estimate(data_train, options);

% perform prediction
data_pred = D.predict(theta, data_test);
label = data_test(end,:)
data_pred
disp(['percentage of correct classification : ' ...
    num2str(sum((label-data_pred)==0) * 100 /30)]);