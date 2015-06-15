function test_mxe_mergeinfo

info = struct('iter',0,'a',1);
last = 1;
newinfo = struct('b', {2 3 4}, 'iter', {0 1 2});
maxiter = 10;
[newinfo, newlast] = mxe_mergeinfo(info, last, newinfo, maxiter, true);