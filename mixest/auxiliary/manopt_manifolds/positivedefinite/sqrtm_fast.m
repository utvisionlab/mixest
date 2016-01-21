%% |sqrtm_fast|
% Calculate matrix square root using Schur decomposition
%
% *Syntax*
%
%   As = sqrtm_fast(A)
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Original author: Reshad Hosseini
%
% Change log: 
%

function As=sqrtm_fast(A)

%  ToDO:For Hermitian positive definite there is a faster version based on
%      Chapter 6 of "Functions of Matrices" by N. J. Higham
%
if any(isinf(A(:))) || any(isnan(A(:))) 
    As = eye(size(A));
    return
end
    
[U,T]=schur(A,'complex');
if isequal(T,diag(diag(T))) 
    % case of symmetric or hermitian input
    Ts = diag(sqrt(diag(T)));   
elseif isreal(T) &&  all(diag(T)>0)
    % case of non-symmetric or non-hermitial positive definite input
    Ts = sqrtm_triu_real(T);
elseif isreal(T)
    % case of real input
    Ts = sqrtm_triu_negative(T);
else
    Ts = sqrtm_triu_complex(T);
end
As = U * Ts * U';
if isreal(A)
    if all(imag(As) == 0)
        As = real(As);
    end
end