function testmle(c,k,n)
% TESTMLE

p = size(c,1);
x=vsamp(c, k, n);
y=sum(x,1);
tmp = norm(y);
rbar = tmp/n;
mysolve(p, rbar)
y = y/tmp;
y
dot(y,c)
