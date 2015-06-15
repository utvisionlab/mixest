function samples = vsamp(center, kappa, n)
% VSAMP returns a set of points sampled from a vMF
% 
% SAMPLES = VSAMP(CENTER, KAPPA, N)   returns a sample of N points
% from a multi-dimensional von Mises Fisher distribution having
% mean direction given by CENTER and concentration parameter given
% by KAPPA.
% 
% CENTER is a column vector
% This implementation is based on the algorithm VM* of A.T.A. Wood
% A. T. A. Wood. Simulation of the von-Mises Fisher
% Distribution. Comm Stat. Comput. Simul. (1994) v. 23 157-164
%
% NOTE: Current Limitation disallows center(1,1) = 0
% 
% See also BETARND, MLE

% d > 1 of course
d = size(center,1);			% Dimensionality
l = kappa;				% shorthand
t1 = sqrt(4*l*l + (d-1)*(d-1));
b = (-2*l + t1 )/(d-1);

x0 = (1-b)/(1+b);

samples = zeros(n,d);

m = (d-1)/2;
c = l*x0 + (d-1)*log(1-x0*x0);

for i=1:n
  t = -1000; u = 1;

  while (t < log(u))
    z = betarnd(m , m);			% z is a beta rand var
    u = rand(1);				% u is unif rand var
    w = (1 - (1+b)*z)/(1 - (1-b)*z);
    t = l*w + (d-1)*log(1-x0*w) - c;
  end
  v = unitrand(d-1);			% d-1 dim unit vector from
                                        % unif distr on sphere
  v = v ./ norm(v);
  samples(i,1:d-1) = sqrt(1-w*w)*v';
  samples(i,d) = w;
end

% Now we form A so that we can get Q
% Now we want an orthogonal transf Q such that
% Q*center = [0 0 ... 1]

% A more eff. way is using householder transformations
% Added on: 1/28/03
% A = eye(d);
% A(:,1)=center;
% Q=mgs(A);
% tmp = Q(:,d);				
% Q(:,d) = Q(:,1);
% Q(:,1) = tmp;

[v,b]=house(center);
Q = eye(d)-b*(v*v.');

for i=1:n
  tmpv = Q*samples(i,:)';
  samples(i,:) = tmpv';
end