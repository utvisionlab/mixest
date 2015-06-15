%% |lbfgs_update|
% *Note:* This is a private function.
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Original author: Reshad Hosseini
%
% Change log: 
%

function  [grad_diff_all, desc_dir_all, x_all, gd_all, ddgd_all, Hdiag, ...
    Expc_all, Expci_all] ...
    = lbfgs_update(x, M, grad_diff, desc_dir, grad_diff_all, ...
    desc_dir_all,num_corr, x_all, gd_all, ddgd_all, Hdiag, ...
    Expc, Expci, Expc_all, Expci_all)
%
% function  [grad_diff_all, desc_dir_all, x_all, gd, ddgd, Hdiag] = ...
%    mbfgs_update(x, M, grad_diff, desc_dir, grad_diff_all, ...
%    desc_dir_all,num_corr, x_all, gd, ddgd, Hdiag)
%
% This function implements Riemmanian BFGS update algorithm
%
% Inputs:
%             x: new point in the manifold
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
%

gd = M.inner(x,desc_dir,grad_diff);
if gd > 1e-10
    num = length(desc_dir_all);
    dd = M.inner(x,desc_dir, desc_dir);
    gg = M.inner(x,grad_diff, grad_diff);
    ddgd = dd / gd;
    if num < num_corr
        desc_dir_all{num+1} = desc_dir;
        grad_diff_all{num+1} = grad_diff;
        x_all{num+2} = x;
        gd_all{num+1} = gd;
        ddgd_all{num+1} = ddgd;
        if isfield(M,'transpf')
            Expc_all{num+1} = Expc;
            Expci_all{num+1} = Expci;
        end
    else
        desc_dir_all = {desc_dir_all{2:num}  desc_dir};
        grad_diff_all = {grad_diff_all{2:num}  grad_diff};
        x_all = {x_all{2:num+1} x};
        gd_all = {gd_all{2:num} gd};
        ddgd_all = {ddgd_all{2:num} ddgd};
        if isfield(M,'transpf')
            Expc_all = {Expc_all{2:num} Expc};
            Expci_all = {Expci_all{2:num} Expci};
        end
    end

    % Initial Hessian
    Hdiag = gd / gg;
else
    if 1
        warning('inner(gradient,descend) is small ... Remove Memory\n');
    end
    % Cleaning History
    desc_dir_all = {};
    grad_diff_all = {};
    x_all = {x};
    gd_all = {};
    ddgd_all = {};
    if isfield(M,'transpf')
        Expc_all = {};
        Expci_all = {};
    end
    %Hdiag = M.norm(x, desc_dir);
    Hdiag = 1;
end