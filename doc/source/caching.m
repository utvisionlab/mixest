%% Caching Intermediate Results

%%
% Some distribution functions allow for a special input and output named
% |store|. This is a structure wherein the function can store the
% intermediate calculation results that it has done for a set of inputs, as
% its fields. These results might speed up things when calling another
% function of the distribution with the same set of inputs. For instance,
% the usual syntax for the |ll| function which calculates log-likelihoods
% is:
%
%   ll = D.ll(theta, data)
%
% But this function also admits the following syntax:
%
%   [ll, store] = D.ll(theta, data, store)
%
% where, either |store| is optional. If we also need the gradient of the
% log-likelihood (|llgrad|) for the same |theta| and |data|, we can use the
% |store| output from the |ll| function as an input to the |llgrad|
% function to prevent |llgrad| from re-calculating the intermediate values
% which are shared in the two functions:
%
%   [ll, store] = D.ll(theta, data);
%   dll         = D.llgrad(theta, data, store);
%
% Please note that you should only share the |store| while working with the
% same distribution and inputs (|theta| and |data|). When these are
% changed, a new store structure should be used, or you will get wrong
% results.
%
