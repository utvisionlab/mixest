%% |mxe_addsharedfields|
% *Note:* This is a private function.
%
% Add shared fields to distribution structures
%
% *Syntax*
%
%   D = mxe_addsharedfields(D)
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Contributors:
%  Reshad Hosseini
%  Mohamadreza Mash'al
%
% Change log: 
%

function D = mxe_addsharedfields(D, notshared)
% Note: fields added here should be documented in doc_distribution_common.m
% and referenced in every distribution factory

    if nargin < 2
        notshared = struct;
    end

% |estimate|
    D.estimate = @estimate;
    function [varargout] = estimate(varargin)
        [varargout{1:nargout}] = mxe_estimate(D, varargin{:}); 
    end

% |randparam|
    D.randparam = @() D.M.rand(); % generate random parameters

% basic operations on Riemannian gradients
    D.sumgrad = @sumgrad;
    function rgrad = sumgrad(rgrad1, rgrad2, theta)
        if nargin < 3
            theta = [];
        end
        rgrad = D.M.lincomb(theta, 1, rgrad1, 1, rgrad2);
    end

    D.scalegrad = @scalegrad;
    function rgrad = scalegrad(scalar, rgrad, theta)
        if nargin < 3
            theta = [];
        end
        rgrad = D.M.lincomb(theta, scalar, rgrad);
    end

% Information Criteria
    D.AICc = @(data) mxe_AICc(D, data);
    D.BIC = @(data) mxe_BIC(D, data);
    
% |pdf|
    if ~isfield(notshared, 'pdf')
        D.pdf = @pdf;
    end
    function y = pdf(theta, data)
        y = exp( D.llvec(theta, data) );
    end
end
