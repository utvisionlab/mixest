function [data,labels] = emsamp(mixture, n)
% EMSAMP    Generates a sample of N data points from Movmf in mixture
% 
% Description
% DATA = EMSAMP(MIXTURE,N) generates a data matrix that has N rows
% (n data points) and each point is drawn from the Movmf MIXTURE 
% [DATA, LABELS] = EMSAMP(MIXTURE,N) generates along with the data
% a labeling from which gaussian in the mixture a point was
% generated. 
% Determine number to sample from each component. Now i am not sure
% if this function is wokring properly?????
% 
% SEE also VSAMP
% hacked over from i think from the NN toolbox..
%
% Author: Suvrit Sra
% (c) The University of Texas at Austin


labels=zeros(n,1);
priors = rand(1,n);

% Pre-allocate data array
data = zeros(n, mixture.dim);

cum_prior = 0;		% Cumulative sum of priors
total_samples = 0;	% Cumulative sum of number of sampled points
for j = 1:mixture.num_clus
  num_samples = sum(priors >= cum_prior & ...
    priors < cum_prior + mixture.priors(j));
  kappa = mixture.kappas(j);
  data(total_samples+1:total_samples+num_samples, :) = ...
      vsamp((mixture.centers(j, :))', kappa, num_samples);
  labels(total_samples+1:total_samples+num_samples) = j;
  cum_prior = cum_prior + mixture.priors(j);
  total_samples = total_samples + num_samples;
end
