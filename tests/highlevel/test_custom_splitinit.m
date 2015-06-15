function test_custom_splitinit

    clc
    clear
    close all

    load data2d

    D = mixturefactory(mvnfactory(2), 1);
    theta = D.init(data);

    options.sm.splitinit = @splitinit;
    
    D.visualize(theta, struct('data', data));

    [newD, newtheta] = D.split(1, theta, options, data);

    newD.visualize(newtheta, struct('data', data));

end

function [newtheta, store] = splitinit(D, idx, theta, options, data, store)
    
    % add parameters for an additional component
    newtheta = theta;
    newtheta.D{end+1} = theta.D{idx};
    % initialize the weights
    newtheta.p(idx) = theta.p(idx) / 2;
    newtheta.p(end+1) = newtheta.p(idx);

    % initialize the means
    Didx = D.component(idx);
    pp = Didx.llvec(theta.D{idx}, data);
    [~, I] = sort(pp, 'descend');
    I = I(1:10);
    [~, C] = kmeans(data(:,I).', 2);
    newtheta.D{idx}.mu = C(1,:).';
    newtheta.D{end}.mu = C(2,:).';
    
    % initialize the covariance matrices
    A = theta.D{idx}.sigma;
    d = size(A,1);
    newtheta.D{idx}.sigma = det(A)^(1/d) * eye(d);
    newtheta.D{end}.sigma = newtheta.D{idx}.sigma;
end
