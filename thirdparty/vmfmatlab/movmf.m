function [clust] = movmf(vectors,k,varargin)
% MOVMF -- Clusters using mixture of vMF distributions
%
% ---------------------------------------------------------
% CLUST = MOVMF(VECTORS, K)
%    VECTORS matrix of data points, each row is a datapoint
%            these are the vectors to be clustered
%    K       number of clusters desired
%
% CLUST = MOVMF(VECTORS, K, TRUTH)
%    VECTORS the input data points
%    K       the number of clusters
%    TRUTH   true cluster label of each data point
%            each entry is in {1,...,k}
% 
% ---------------------------------------------------------
% Author: Arindam Banerjee
% Minor modifications by: Suvrit Sra 
%
% Copyright lies with the author
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
%
% ---------------------------------------------------------

if (nargin == 3)
  truth=varargin{1};
end

[D,V] = size(vectors);
dim   = V;

% --------------------------------
% Getting the initial random means
% --------------------------------
% getting global mean
sumv  = sum(vectors);
normv = sqrt(sumv*sumv');
mu0   = sumv./normv;

% perturbing global mean to get initial cluster centroids
perturb = 0.1;
for h = 1:k
  randVec     = rand(1,V) - 0.5;
  randNorm    = perturb*rand;

  smallRandVec =  randNorm*randVec/sqrt(randVec*randVec'); 
  mu(h,:)      =  mu0 + smallRandVec;
  mu(h,:)      =  mu(h,:)/sqrt(mu(h,:)*mu(h,:)');
end

% --------------------------------
% Getting the means from spkmeans
% --------------------------------
diff    = 1;
epsilon = 0.001;
value   = 100;
iteration = 1;
while (diff > epsilon)
  
  display(['Iteration ',num2str(iteration)]);
  iteration = iteration+1;
  oldvalue      = value;
  
  % assign points to nearest cluster
  simMat        =  vectors*mu';
  [simax,clust] =  max(simMat,[],2);

  % compute objective function value
  value         = sum(simax);
  
  % compute cluster centroids
  for h=1:k
    sumVec(h,:)  = sum(vectors(find(clust==h),:));
    mu(h,:)      = sumVec(h,:)/sqrt(sumVec(h,:)*sumVec(h,:)');
  end
  
  diff = abs(value - oldvalue);
  
end

%display(clust);
figure;
subplot(2,1,1),plot(1:D,clust,'bo');
display('Initial iterations done');

Clust1 = clust;

%----------------------------------------------
% You can cut the code at this point
%----------------------------------------------

% --------------------------------
% movMF iterations
% --------------------------------

diff      = 1;
epsilon   = 0.0001;
value     = 100;
iteration = 1;

kappaMax = 5000;
kappaMin = 1;

% initializing kappa, alpha
kappa    = kappaMin*ones(1,k);
for h = 1:k
  alpha(h) = length(find(clust==h))/D;
end

T = 100;
t = 1;

display('Starting main iterations ...');
while (diff > epsilon)
%while (iteration < 10)

  display(['Iteration ',num2str(iteration)]);
  iteration = iteration + 1;
  oldvalue   = value;
  
  % assignment of points to nearest vMFs
  
  logNormalize  = log(alpha) + (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(dim/2-1,kappa); 
  logProbMat    = (vectors*(mu'.*(ones(dim,1)*kappa)) + ones(D,1)* ...
      logNormalize)/T; 
  lpm = logProbMat;
  logSum        = log(sum(exp(logProbMat),2)); % this step, without
                                               % the 1/T in the
                                               % previous line, leads to
                                               % Inf. We do it w/o
                                               % 1/T in C++ using NTL
					       
  logProbMat    = logProbMat - logSum*ones(1,k);
  
  value = sum(sum(logProbMat));
  
  
  % updating component parameters
  alpha  = sum(exp(logProbMat)); % this step leads to Inf
  mu     = exp(logProbMat')*vectors;
  
  for h=1:k
    oldkappa(h) = kappa(h);

    normMu   = sqrt(mu(h,:)*mu(h,:)');
    rbar(h)  = normMu/(t*alpha(h));
    
    mu(h,:)  = mu(h,:)/normMu;
    kappa(h) = (rbar(h)*dim - rbar(h)^3)/(1-rbar(h)^2);
    alpha(h) = alpha(h)/D;
    display(['Size = ',num2str(D*alpha(h)),' rbar =',num2str(rbar(h)),' kappa = ',num2str(kappa(h))]);
  end

  diff = abs(value - oldvalue);
  
end


[simax,clust] = max(logProbMat,[],2);  
subplot(2,1,2),plot(1:D,clust,'ro');
Clust2 = clust;


if (nargin == 3)
  % ----
  % evaluation using confusion matrices

  c = length(unique(truth));

  conf1 = zeros(k,c);
  conf2 = zeros(k,c);
  for i=1:D
    conf1(Clust1(i),truth(i)) =  conf1(Clust1(i),truth(i)) + 1;
    conf2(Clust2(i),truth(i)) =  conf2(Clust2(i),truth(i)) + 1;
  end

  display('Results from Initial Convergence');
  display(conf1);

  display('Results from movMF Convergence');
  display(conf2);
end