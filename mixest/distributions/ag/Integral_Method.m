function A = Integral_Method(mat, Tol)
%%
% Calculating the expected logarithm of weighted sum of chi-square
% distribution (A = E[ \log \sum d_i N_i^2]) with the integral method
%
% Written by: Pourya H. Zadeh

% eigs should be a column vector
[r, c] = size(mat);
if r == c
    eigs = eig(mat);
else
    eigs = mat;
end

if size(eigs,2) > 1
    eigs = eigs.';
end

if nargin == 1
    Tol = 1e-13;
end

eigs = sort(eigs);
d = length(eigs);

fun = @(x) inint(x, eigs, d);
cons = quadgk(fun,0,Inf,'RelTol',0,'AbsTol',Tol);
A = psi(d/2)+log(2*sum(eigs)/d)+cons;
end

function y = inint(x, d, n)
%% Function for value inside integral
% x should have one row and others should have one column
im = 1;
for j = 1:n
    im = im .* ((1+2*d(j).*x).^(-1/2));
end

y = (1./x).*(((1+(2*x*sum(d)/n)).^(-n/2))-im);
end