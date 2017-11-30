function M = mxe_addsharedmanifold(M)
% *Note:* This is a private function.
%
% Add shared fields to manifold structures
% This is needed for applying LBFGS method on manifolds
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

% adding store to the last column of retraction
retr = M.retr;
retr2 = @retraction;
    function [varargout, store] = retraction(varargin)
        if nargin > 3
            store = varargin{4};
            varargin = varargin(1:3);
        end
        varargout = retr(varargin{:});
    end
    if nargin(M.retr) == 3
        M.retr = retr2;
    end

% adding store to the last column of vector transport
transp = M.transp;
transp2 = @transpvec;
    function [matout, store] = transpvec(X, Y, E, varargin)
        if nargin > 5
            store = varargin{3};
        end
        %varargin = varargin(1:3);
        matout = transp(X, Y, E);
    end
    if nargin(M.transp) == 3
        M.transp = transp2;
    end
    
% adding operation of retraction and transp together
retrtransp = @retractiontranspvec;
    function [Y, F, store] = retractiontranspvec(X, U, E, t, store)
        if nargin < 5
            store = struct;
        end
        if nargin < 4
            t = 1;
        end
        [Y, store] = M.retr(X, U, t, store);
        [F, store] = M.transp(X, Y, E, U, t, store);
    end
if ~isfield(M, 'retrtransp')
    M.retrtransp = retrtransp;
end

% adding operation of retraction and transp together
transpdiffE = @transpvecdifferentE;
    function [F, store] = transpvecdifferentE(X, U, E, t, store)
        if nargin < 5
            store = struct;
        end
        if nargin < 4
            t = 1;
        end
        [Y, store] = M.retr(X, U, t, store);
        [F, store] = M.transp(X, Y, E, U, t, store);
    end
if ~isfield(M, 'transpdiffE')
    M.transpdiffE = transpdiffE;
end

% adding fast version of transpvec
transpstore = @transpvecstore;
    function [expconstruct, iexpconstruct] = transpvecstore(X, Y)
        expconstruct.X = X;
        expconstruct.Y = Y;
        iexpconstruct = expconstruct;
    end
if ~isfield(M, 'transpstore')
    M.transpstore = transpstore;
end

% faster version of vector transport by storing some information
transpf = @transpvecfast; 
    function F = transpvecfast(expconstruct, E)
        %expconstruct
        F = M.transp(expconstruct.X, expconstruct.Y, E);
    end
if ~isfield(M, 'transpf')
    M.transpf = transpf;
end

% faster version of adjoint vector transport by storing some information
atranspf = @atranspvecfast; 
    function F = atranspvecfast(iexpconstruct, E)
        F = M.transp(iexpconstruct.Y, iexpconstruct.X, E);
    end
if ~isfield(M, 'atranspf')
    M.atranspf = atranspf;
end

end
