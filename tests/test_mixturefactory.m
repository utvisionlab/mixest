function test_suite = test_mixturefactory %#ok<STOUT>
initTestSuite;

function s = setup %#ok<*DEFNU>
% This is run before running each test case and the output is passed to the
% test case.

s.D = mixturefactory(mvnfactory(1), 5);

s.theta.D{1} = struct('mu',1 ,'sigma',1);
s.theta.D{2} = struct('mu',2 ,'sigma',2);
s.theta.D{3} = struct('mu',3 ,'sigma',3);
s.theta.D{4} = struct('mu',4 ,'sigma',4);
s.theta.D{5} = struct('mu',5 ,'sigma',5);
s.theta.p = [0.1; 0.2; 0.3; 0.12; 0.28];

% 10 samples from 0.3*N(-1,1) + 0.7*N(5,4)
s.data = [4.873890253620687,6.429485807652192,-2.349886940156521,4.590067883400451,4.751711303567376,2.034923466331855,-0.274595775053894,7.979395215570930,7.818068979600959,7.834384826859228];




function test_num(s)
assert(s.D.num() == 5)
D = s.D.fixate([4,5], s.theta);
assert(D.num() == 3)


function test_numfixed(s)
assert(s.D.numfixed() == 0)
D = s.D.fixate([4,5], s.theta);
assert(D.numfixed() == 2)


function test_numtotal(s)
assert(s.D.numtotal() == 5)
D = s.D.fixate([4,5], s.theta);
assert(D.numtotal() == 5)


function test_component(s)
D2 = s.D.component(2);
assert(D2.datadim() == 1)
D23 = s.D.component([2,3]);
assert(numel(D23) == 2)


function test_varD(s)
D3 = s.D.varD(3);
assert(D3.datadim() == 1)
D23 = s.D.varD([2,3]);
assert(numel(D23) == 2)


function test_fixedD(s)
D = s.D.fixate([4,5], s.theta);
D2 = D.fixedD(2);
assert(D2.datadim() == 1)
D12 = D.fixedD([1,2]);
assert(numel(D12) == 2)


function test_fixedparam(s)
D = s.D.fixate([4,5], s.theta);
fixedtheta = D.fixedparam();
assertElementsAlmostEqual(fixedtheta.p, [0.3; 0.7])


function test_subparam(s)
theta = s.D.subparam(s.theta, 2:3);
assert(numel(theta.D) == 2)
assert(theta.D{2}.mu == 3)


function test_fullparam(s)
[D, theta] = s.D.fixate([4,5], s.theta);
fulltheta = D.fullparam(theta);
assert(fulltheta.D{2}.mu == 2)
assert(fulltheta.D{5}.mu == 5)
assertElementsAlmostEqual(fulltheta.p, [0.1; 0.2; 0.3; 0.12; 0.28])


function test_addcomponent(s)
newD = s.D.addcomponent(mvnfactory(1));
assert(newD.num() == 6)


function test_removecomponent(s)
newD = s.D.removecomponent(2);
assert(newD.num() == 4)


function test_invertindex(s)
invidx = s.D.invertindex([2,4]);
assertEqual(invidx, [1,3,5])
newD = s.D.fixate([4,5], s.theta);
invidx = newD.invertindex(2, 'fixed');
assertEqual(invidx, 1)


function test_fixate(s)
[D, theta] = s.D.fixate([4,5], s.theta);
[newD, newtheta] = D.fixate(2:3, theta);
assert(D.num() == 3)
assert(newD.num() == 1)
assert(newD.numfixed() == 4)
assertElementsAlmostEqual(newtheta.p, [0.1; 0.9])
fixedTheta = newD.fixedparam();
assertElementsAlmostEqual(fixedTheta.p, [0.133333333333333;0.311111111111111;0.222222222222222;0.333333333333333])

randtheta = newD.randparam();
fulltheta = newD.fullparam(randtheta);
assert(fulltheta.D{5}.mu == 3)

newD = s.D.fixate('all', s.theta);
assert(newD.num() == 0)
assert(newD.numfixed() == 5)



function test_unfix(s)
[D, theta] = s.D.fixate([4,5], s.theta);
[newD, newtheta] = D.unfix(2, [], theta);
assert(D.num() == 3)
assert(newD.num() == 4)
assert(newD.numfixed() == 1)
assert(newtheta.D{4}.mu == 5)
assertElementsAlmostEqual(newtheta.p, [0.1; 0.2; 0.3; 0.28; 0.12])

newD = s.D.unfix('all');
assert(newD.num() == 5)
assert(newD.numfixed() == 0)


function test_split(s)
[D, theta] = s.D.fixate([4,5], s.theta);
[newD, newtheta, idxSplitted, idxMap] = D.split(2, theta);
assert(newD.num() == 4)
assert(newD.numfixed() == 2)
assertElementsAlmostEqual(newtheta.p, [0.1; 0.1; 0.3; 0.1; 0.4])
assert(numel(newD.varD(idxSplitted)) == 2)
idxUnchanged = D.invertindex(2);
assertEqual(idxMap(idxUnchanged), [1,3])


function test_merge(s)
[D, theta] = s.D.fixate([4,5], s.theta);
[newD, newtheta, idxMerged, idxMap] = D.merge(1, 2, theta);
assert(newD.num() == 2)
assert(newD.numfixed() == 2)
assertElementsAlmostEqual(newtheta.p, [0.3; 0.3; 0.4])
assert(numel(newD.varD(idxMerged)) == 1)
idxUnchanged = D.invertindex([1,2]);
assertEqual(idxMap(idxUnchanged), 2)


function test_dim(s)
assert(s.D.dim() == 15) % 5(num) * 2(D) + 5(p)
D = s.D.fixate([4,5], s.theta);
assert(D.dim() == 10) % 3(num) * 2(D) + 4(p)


function test_datadim(s)
assert(s.D.datadim() == 1)
D = s.D.fixate([4,5], s.theta);
assert(D.datadim() == 1)


function test_ll(s)
ll = s.D.ll(s.theta, s.data);


function test_llvec(s)
llvec = s.D.llvec(s.theta, s.data);


function test_llgrad(s)
dll = s.D.llgrad(s.theta, s.data);


function test_llgraddata(s)
dld = s.D.llgraddata(s.theta, s.data);


function test_cdf(s)
y = s.D.cdf(s.theta, s.data);


function test_pdf(s)
y = s.D.pdf(s.theta, s.data);


function test_sample(s)
data = s.D.sample(s.theta);
assertEqual(size(data), [1,1]);

data = s.D.sample(s.theta, 5);
assertEqual(size(data), [1,5]);


function test_entropy(s)
h = s.D.entropy(s.theta);


function test_kl(s)
kl = s.D.kl(s.theta);


%TODO
% function test_penalizer(s)
% [costP, gradP] = s.D.penalizer(s.theta, s.data);
% assertElementsAlmostEqual(costP, 0);
% assertElementsAlmostEqual(gradP.mu, 0);
% assertElementsAlmostEqual(gradP.sigma, 0);


% function test_estimatedefault(s)
% options.crossval = true;
% theta = s.D.estimatedefault(s.data, options);
% 
% 
% function test_estimatepartial(s)
% options.crossval = true;
% theta = s.D.estimatepartial(1:2, s.theta, s.data, options);


function test_init(s)
theta = s.D.init(s.data);


function test_randparam(s)
theta = s.D.randparam();


function test_sumparam(s)
theta = s.D.sumparam(s.theta, s.theta);


function test_scaleparam(s)
theta = s.D.scaleparam(-1, s.theta);


function test_sumgrad(s)
grad1 = s.D.M.egrad2rgrad(s.theta, s.theta);
grad = s.D.sumgrad(grad1, grad1, s.theta);


function test_scalegrad(s)
grad1 = s.D.M.egrad2rgrad(s.theta, s.theta);
grad = s.D.scalegrad(-1, grad1, s.theta);


%TODO test information criteria functions
