---
title: 'MixEst: An Estimation Toolbox for Mixture Models'
tags:
  - mixture models 
  - mixtures of experts 
  - manifold optimization
  - expectation-maximization 
  - stochastic optimization
authors:
 - name: Reshad Hosseini
   orcid: 0000-0002-3669-760X
   affiliation: 1
 - name: Mohamadreza Mash'al
   affiliation: 1
affiliations:
 - name: School of ECE, College of Engineering, University of Tehran
   index: 1
date: 15 December 2017
bibliography: paper.bib
---

# Summary #
    
Mixture models are powerful statistical models used in many applications ranging from density estimation to clustering and classification.
When dealing with mixture models, there are many issues that the experimenter should be aware of and needs to solve.
The MixEst toolbox is a powerful and user-friendly package for MATLAB that implements several state-of-the-art approaches to address these problems.
Additionally, MixEst gives the possibility of using manifold optimization for fitting the density model, a feature specific to this toolbox.
MixEst simplifies using and integration of mixture models in statistical models and applications.
For developing mixture models of new densities, the user just needs to provide a few functions for that statistical distribution and the toolbox takes care of all the issues regarding mixture models.

# Introduction #

Mixture models are an integrated and fundamental component in many
machine learning problems ranging from clustering to regression and
classification [@McLPee00]. Estimating the parameters of mixture models
is a challenging task due to the need to solve the following issues in
mixture modeling:

-   Unboundedness of the likelihood: This problem occurs when one
    component gets a small number of data points and its likelihood
    becomes infinite [@ciuperca_penalized_2003].

-   Local maxima: The log-likelihood objective function for estimating
    the parameters of mixture models is non-concave and has many local
    maxima [@ueda_split_2000].

-   Correct number of components: In many applications, it is needed to
    find the correct number of components [@khalili_variable_2007].

Addressing these issues for a mixture density when it is not available
in common mixture modeling toolboxes will cost a lot of time and effort
for the experimenter. MixEst addresses all these issues not only for
already implemented densities, but also for densities that the user may
implement. By implementing densities, we mean implementing a few simple
functions which will be briefly discussed in [Model Development Section](#model-development).

This toolbox provides a framework for applying manifold optimization for
estimating the parameters of mixture models. This is an important
feature of this toolbox, because recent empirical evidence shows that
manifold optimization can surpass expectation maximization in the case
of mixtures of Gaussians [@hosseini2015]. It also opens the door for
large-scale optimization by using stochastic optimization methods.
Stochastic optimization also allows solving the likelihood unboundedness
problem mentioned above, without the need of implementing a penalizing
function for the parameters of the density.

While several libraries are available for working with mixture models,
to the best of our knowledge, none of them offers a modular and flexible
framework that allows for fine-tuning the model structure or can provide
universal algorithms for estimating model parameters solving all the
problems listed above. A review of features available in some libraries
can be seen in [Comparison
Section](#feature-comparison).

In the next section, we give a short overview of the toolbox and its
features.

# About the MixEst Toolbox # 

This toolbox offers methods for constructing and estimating mixtures for
joint density and conditional density modeling, therefore it is
applicable to a wide variety of applications like clustering, regression
and classification through probabilistic model-based approach. Each
distribution in this toolbox is a structure containing a manifold
structure representing parameter space of the distribution along with
several function handles implementing density-specific functions like
log-likelihood, sampling, etc. Distribution structures are constructed
by calling factory functions with some appropriate input arguments
defining the distribution. For example for constructing a mixture of
one-dimensional Gaussians with 2 components, it will suffice to write
the following commands in MATLAB:

    Dmvn = mvnfactory(1);
    Dmix = mixturefactory(Dmvn, 2);

As an example of how to evoke a function handle, consider generating
1000 samples from the previously defined mixture:

    theta.D{1}.mu = 0; theta.D{1}.sigma = 1; % mean and variance of the 1st component
    theta.D{2}.mu = 5; theta.D{2}.sigma = 2; % mean and variance of the 2nd component
    theta.p = [0.8 0.2]; % weighting coefficients of components
    data = Dmix.sample(theta, 1000);

Each distribution structure exposes a common interface that optimization
algorithms in the toolbox can use to estimate its parameters. In
addition to the EM algorithm which is a commonly implemented method in
available libraries, our toolbox also makes optimization on manifolds
available featuring procedures like early-stopping and mini-batching to
avoid overfitting. For optimization on manifolds, our toolbox depends on
optimization procedures of an excellent toolbox called Manopt
[@boumal_manopt_2014]. In addition to optimization algorithms of Manopt
like steepest descent, conjugate gradient and trust regions methods, the
user can also use our implementation of Riemmanian LBFGS method.

# Model Development #

MixEst includes many joint and conditional distributions to model data
ranging from continuous to discrete and also directional. Some users,
however, may want to apply the tools developed in this toolbox for
mixtures of a distribution not available in the toolbox yet. To this
end, the user needs to write a factory function that constructs a
structure for the new distribution.

Each distribution structure has a field named "M\" determining the
manifold of its parameter space. For example for the case of
multivariate Gaussian distribution, this is a product manifold of a
positive definite manifold and a Euclidean manifold:

    % datadim is the function input argument determining the dimensionality of data
    muM = euclideanfactory(datadim);
    sigmaM = spdfactory(datadim);
    D.M = productmanifold(struct('mu', muM, 'sigma', sigmaM));

The manifold of parameter space completely determines how parameter
structure is given to or is returned by different functions. The
structure of parameters for multivariate Gaussian would have two fields,
a mean vector "mu\" and a covariance matrix "sigma\".

To use the estimation tools of the toolbox, two main functions have to
be implemented. The *weighted log-likelihood* (wll) function and a
function for computing the gradient of sum-wll with respect to the
distribution parameters. The syntax for calling the wll function is:

    llvec = D.llvec(theta, data);

The input argument `theta` is a structure containing the input
parameters of the corresponding distribution. The second input argument
`data` can be either a data matrix or a structure having several fields
such as the data matrix and weights, which is interpreted using the
`mxe_readdata` function. The output argument `llvec` is a vector with
entries equal to wll for each datum (each column) in the data matrix.

The function to compute the gradient of sum-wll has the following
syntax:

    llgrad = D.llgrad(theta, data);

The input arguments are similar to the function `llvec`. The output
argument `llgrad` is a structure similar to the input argument `theta`
returning the gradient of sum-wll with respect to each parameter.

Some other (optional) functions that can be implemented for
distributions are:

-   `init`: This is for initializing the estimator using the data.

-   `estimatedefault`: If the maximum wll has a structure that allows
    fast optimization (or has a closed-form solution), this estimator
    can be implemented in this function. When this function is not
    present, the Riemmanian optimization is called in the maximization
    step of EM algorithm.

-   `llgraddata`: This function computes the gradient of wll with
    respect to the data. It is required in some special cases such as
    when the distribution is used as the radial component of an
    elliptically-contoured distribution or as the components in
    independent component analysis.

-   `ll`: This function is sum-wll (sum of the output vector of `llvec`
    function). Sometimes it is faster to write this function differently
    than just calling `llvec` and summing up its output vector.

Two other functions that can be used in the split-and-merge algorithms
to avoid local maxima of mixture models are `kl` (for computing
KL-divergence) and `entropy` (for computing entropy). If the user wants
to evoke a maximum-a-posteriori estimate, the functions
`penalizerparam`, `penalizercost` and `penalizergrad` need to be
implemented.

# Feature Comparison #

To demonstrate the richness of features in MixEst, we are comparing its
features with several other well-known packages in
The following table. Among many toolboxes available for mixture
modeling, we select those that are feature-rich and representative.
These packages are Sklearn [@scikit-learn], Mclust [@fraley1999mclust],
FlexMix [@Leisch_2004b], Bayes Net [@Murphy01thebayes] and
MixMod [@biernacki2006model]. We include Bayes Net to demonstrate what a
generic Bayesian graphical modeling toolbox can do. Sklearn is a
powerful machine learning toolbox containing many tools, among others
tools specific for mixture modeling. MixMod also provides bindings for
Scilab and Matlab.

         |  MixEst |  SKlearn |  Mclust |  FlexMix |  Bayes Net |  MixMod
  -------| --------| ---------| --------| ---------| -----------| --------
  Programming language   |  Matlab |  Python  |    R    |     R    |   Matlab   |   C++
  Approaches for solving local minima   |    SM   |   IDMM   |    HC   |    ---   |     ---    |   ---
  Manifold optimization   |   Yes   |    No    |    No   |    No    |     No     |    No
  Bayesian approaches for inference   |   MAP   |    VB    |   MAP   |    ---   |     MAP    |    SM
  Large-scale optimization   |    MB   |    ---   |   ---   |    ---   |     ---    |   SEM
  Having tools for model selection   |   Yes   |    No    |   Yes   |    No    |     No     |   Yes
  Automatic model selection   |   CSM   |   IDMM   |   ---   |    ---   |     ---    |   ---
  Ease of extensibility   |   Easy  |    ---   |   ---   |   Easy   |   Medium   |   ---
  Having mixtures of experts   |   Yes   |    No    |    No   |    No    |     Yes    |    No
  Having mixtures of classifiers  |   Yes   |    No    |    No   |    No    |     Yes    |    No
  Having mixtures of regressors  |   Yes   |    No    |    No   |    Yes   |     Yes    |    No

[Table:  SM stands for split-and-merge approach,
  IDMM stands for infinite dirichlet mixture models, HC stands for
  initialization using hierarchical clustering; MAP stands for
  maximum-a-posteriori, VB stands for variational Bayes; SEM stands for stochastic EM, MB stands for  mini-batching; CSM stands for competitive split-and-merge]
  
# References ##
