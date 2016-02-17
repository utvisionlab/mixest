function [X,T,L1,L2] = mixspheregen(n,m,k,d,c,e)

%mixgen - Gaussian mixture generator 
%
%[X,T] = mixgen(n,m,k,d,c,e) 
%  n - size of training set 
%  m - size of test set 
%  k - number of components
%  d - dimension
%  c - separation degree
%  e - maximum ratio of maximum kappa to minimum kappa
%returns
%  X - training set (n x d)
%  T - test set (m x d) 

% Reshad Hosseini and Mohammad-Reza Mash'al, 2015
% Modified from Nikos Vlassis, 2000


R=zeros(k,d^2);

% mixing weights
while 1
  W = rand(k,1); 
  W = W / sum(W);
  if all(W > 1/(4*k)) 
    break;
  end
end


stdMax = 0.1;
stdMin = 0.1 / e;

stdvar = ['stdev' num2str(d)];
kappavar = ['kappa' num2str(d)];
if exist('stdkappa','file')
    load('stdkappa', stdvar, kappavar);
end

if ~exist(stdvar,'var') && ~exist(kappavar,'var')
    kappa = logspace(1+log10(d),4.5+log10(d),100);
    stdev = zeros(1,length(kappa));
    for k1=1:length(kappa)
        center = zeros(d, 1);
        center(end) = 1;
        samples = vsamp(center, kappa(k1), 10000);
        dist = acos( samples * center);
        fprintf('.');
        stdev(k1) = sqrt(sum(dist.^2)/ 10000);
    end
    fprintf('\n');
    eval([ stdvar '= stdev']);
    eval([ kappavar '= kappa']); 
    if exist('stdkappa','file')
        save('stdkappa', stdvar, kappavar,'-append');
    else
        save('stdkappa', stdvar, kappavar);
    end
end

eval([ 'stdev = ' stdvar ';']);
eval([ 'kappa = ' kappavar ';']); 

% computing kappa from std
std = stdMin + rand(k,1) * (stdMax - stdMin);
kap = interp1(stdev, kappa , std);

% create c-separated vMF clusters of maximum ratio e
trials = 0;

x = zeros(d, k);
x(end, :) = 1;
while 1 
  rd = randn(d, k);
  rd(end, :) = 0;
  nd = sqrt(sum(rd.^2, 1));
  rd = bsxfun(@rdivide, rd, nd);
  
  t = sqrt(c) * sqrt(k) * trials/10 * rand(1,k);
  
  M = bsxfun(@times, x, cos(t)) + bsxfun(@times, rd, sin(t));
  
  % check degree of separation
  error = 0;
  for i = 1:k-1
    for j = i+1:k
      if acos(M(:,i).'*M(:,j)) < c * max(std(i),std(j))
        error = 1;
      end
    end
  end
  if ~error
    break;
  end
  trials = trials + 0.1;
end


mixtures.priors= W;
mixtures.centers = M.';
mixtures.kappas= kap;
mixtures.dim = d;
mixtures.num_clus = k;
fprintf('|.');
X = emsamp(mixtures, n);
fprintf('.\n');
T = emsamp(mixtures, m);
L1 = Inf;
L2 = Inf;