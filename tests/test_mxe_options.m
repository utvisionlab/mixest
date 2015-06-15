function test_mxe_options

clear o
o.solver = 'test1';
o.maxiter = 10;
o.inner.maxiter = 5;

oo = mxe_options(o);
assert(oo.maxiter == 10)
assert(oo.inner.maxiter == 5)

ooo = mxe_options(oo); % this should return immediately
ooo = mxe_options(oo, struct('a',1)); % this should go into work


io = mxe_inneroptions(oo, [], 'aaa');
assertEqual(io.solver, 'test1')

oo.inner.aaa.solver = 'test2';
io = mxe_inneroptions(oo, [], 'aaa');
assertEqual(io.solver, 'test2')

