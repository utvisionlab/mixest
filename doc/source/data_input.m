%% Data Input Argument to Functions
% In MixEst, all the functions taking a |data| input argument, can accept
% the data in one of the following forms:
%

%% Data Matrix
% In matrix form, the data is arranged in an n-by-N matrix where |n| is
% the dimensions of the data space and |N| is the number of data points
% (i.e. each column of the matrix represents a single observation).
%
% *Example*
%
%   % create a distribution and random parameters
%   D = mvnfactory(2);
%   theta = D.randparam();
%   % generate 1000 random 2-dimensional data points
%   data = randn(2, 1000);
%   % pass the data to the ll function
%   ll = D.ll(theta, data)
%
% *Note:* In MATLAB versions R2011b and above, instead of the data matrix,
% you can give a |matlab.io.MatFile| object referring to a MAT-file
% containing a variable named |data| with specifications as pointed out
% above.
%

%% Data Structure
% Using the structure form, you can pass along with the data matrix, the
% indices of a subset of data points and/or a weighting vector on the data
% points whenever applicable. In this form, |data| is a structure with the
% following fields (All fields except |data| are optional):
% 
% * *|data|* (|n-by-N|) The data matrix as described in the previous
% section.
% * *|index|* (|1-by-L|) Data indexing vector for selecting a subset of
% data (Only numeric indexing is supported. Logical indexing is not
% supported).
% * *|weight|* (|1-by-N|, or |1-by-L| if indexing is used) Data weighting
% vector containing values between zero and one to be multiplied by each
% data point.
%
% where |n| is the dimensions of the data space, |N| is the total number of
% data points, and |L| is the number of points in the indexed subset of
% data.
%
% *Example 1*
%
%   % create a distribution and random parameters
%   D = mvnfactory(2);
%   theta = D.randparam();
%   % generate 1000 random 2-dimensional data points
%   data = randn(2, 1000);
%   % use only the 100th up to 200th points from the data
%   index = 100:200;
%   % generate random weights for the selected data
%   weight = rand(1, 101); % 101 is the number of elements in index
%   % build the data structure
%   data = struct('data', data, 'index', index, 'weight', weight);
%   % pass the data to the ll function
%   ll = D.ll(theta, data)
%
% *Example 2*
%
%   D = mvnfactory(2);
%   theta = D.randparam();
%   ll = D.ll(theta, struct('data', randn(2,1000), 'weight', rand(1,1000)))
%