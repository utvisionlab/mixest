%% |desc_dir_cal|
% *Note:* This is a private function.
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Original author: Reshad Hosseini
%
% Change log:
%

function  p_desc = desc_dir_cal(p, M, grad_diff_all, desc_dir_all, ...
    x_all, ddgd, gd, corrections , Hdiag, Expc_all, Expci_all)
%
%function  p_desc = desc_dir_cal(p, M, grad_diff_all, desc_dir_all, ...
%    x_all, gdivg, s_gdivg, s_divg, hy_divg)
%
% This function implements Riemmanian BFGS algorithm for finding the
% descend direction p_desc=H*p , where p is the input vector in the
% tangent space and H is estimated inverse hessian using BFGS algorithm
%
% Inputs:
%             p: A vector in the gradient space of the current point
%
%             M: Manifold object
%
% grad_diff_all: Gradient difference of the previous points
%    grad_diff_all(:,k) = Grad(f,x_{k+1}) - transp(x_k,x_{k+1},Grad(f,x_k))
%
%  desc_dir_all: Descent direction of the previous points
%    desc_dir_all(:,k)= transp(x_k,x_{k+1},alpha*p_desc_prev)
%
%        x_all: All previous points
%
%         ddgd: <desc_dir_k,desc_dir_k>/<grad_diff_k,desc_dir_k>
%
%           gd: <grad_diff_k,desc_dir_k>
%
%  corrections: number of corrections till the point
%
%        Hdiag: The scalar of the initial Hessian
%
% Output: p_desc is the descnedt direction
%

coef = M.inner(x_all{corrections+1}, p, ...
    desc_dir_all{corrections}) / gd{corrections};
p_prev = M.lincomb(x_all{corrections+1}, 1, p , - coef, ...
    grad_diff_all{corrections});

% Doing an inverse retraction from a point to the previous one
if isfield(M,'transpf')
    p_invtransp = M.transpf(Expci_all{corrections},p_prev);
    if corrections >1
        vec_prec= desc_dir_cal(p_invtransp, M, grad_diff_all, desc_dir_all, ...
            x_all, ddgd, gd, corrections-1, Hdiag, Expc_all, Expci_all);
    else
        vec_prec = M.lincomb(x_all{corrections+1}, Hdiag, p_prev);
    end
    vec_new = M.atranspf(Expc_all{corrections},vec_prec);
else
    if isfield(M,'atranspf')
        p_invtransp = M.atranspf(x_all{corrections},x_all{corrections+1},p_prev);
    else
        % inverse vector transport is not implemented
        p_invtransp = M.transp(x_all{corrections+1},x_all{corrections},p_prev);
    end
    if corrections >1
        vec_prec= desc_dir_cal(p_invtransp, M, grad_diff_all, desc_dir_all, ...
            x_all, ddgd, gd, corrections-1, Hdiag);
    else
        vec_prec = M.lincomb(x_all{corrections+1}, Hdiag, p_prev);
    end
    vec_new = M.transp(x_all{corrections},x_all{corrections+1},vec_prec);
end

coef = M.inner(x_all{corrections+1},vec_new,grad_diff_all{corrections});
p_desc = M.lincomb(x_all{corrections+1}, 1, vec_new, ...
    -coef/gd{corrections}, desc_dir_all{corrections});
p_desc =  M.lincomb(x_all{corrections+1}, 1, p_desc, ...
    ddgd{corrections}, p);