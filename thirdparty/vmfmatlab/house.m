function [v,b]=house(x)
% HOUSE Returns the householder transf to reduce x to b*e_n 
%
% [V,B] = HOUSE(X)  Returns vector v and multiplier b so that
% H = eye(n)-b*v*v' is the householder matrix that will transform
% Hx ==> [0 0 0 ... ||x||], where  is a constant.

n=length(x);

s = x(1:n-1)'*x(1:n-1);
v= [x(1:n-1)' 1]';

if (s == 0)
  b = 0;
else
  m = sqrt(x(n)*x(n) + s);
  
  if (x(n) <= 0)
    v(n) = x(n)-m;
  else
    v(n) = -s/(x(n)+m);
  end
  b = 2*v(n)*v(n)/(s + v(n)*v(n));
  v = v/v(n);
end
