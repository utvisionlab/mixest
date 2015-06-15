function sample = unitrand(dim)
% UNITRAND  Generate unit random vector in dim dimensions
% 
% SAMPLE = UNITRAND(DIM)   Generates a unit random vector that
% comes from a uniform distr on the unit 'dim' dimensional sphere.
% The method used is: generate N normally distr numbers, form a
% vector and normalize to be of unit length!
% 
% See: HAKMEN Item 27

sample = randn(dim,1);
n = norm(sample);
sample = sample/n;
