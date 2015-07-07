## About ##

MixEst: A MATLAB toolbox for mixture-model parameter estimation

This toolbox is copyright (C) 2015 by Reshad Hosseini & Mohamadreza Mash'al and is distributed under the terms of the GNU General Public License (GPL) version 3 (or later).

Contact: [Reshad Hosseini](mailto:reshad.hosseini@ut.ac.ir) or [Mohamadreza Mash`al](mailto:mrmashal@ut.ac.ir)

## Quick installation guide ##

* Unzip andd copy the whole MixEst directory you just downloaded in a location of your choice on disk, say, in /my/directory/.

* Go to /my/directory/eg/ at the Matlab command prompt and execute 'install_mixest'. You may save this path for your next Matlab sessions: follow the menu File . Set Path... and save.

##Directory structure##

<pre>
./ The top directory with README.md
 | tests/           - test functions 
 | doc/             - toolbox documentation
 | examples/        - some examples showing how to use MixEst
 | mixest/          - folder containing main toolbox files
   |--- auxiliary/        - functions needed for creating a distribution or estimating its parameters
   |--- distributions/    - a collection of functions for creating distribution structures
   |--- estimation/       - main function for estimating the parameters
   |--- gates/            - gate factories needed for creating mixtures of experts
 | thirdParty/      - third party tools
   |--- lightspeed/       - logsumexp function of T. Minka's matlab toolbox
   |--- randraw/          - function for sampling from several distributions by A. Bar-Guy
   |--- vmfmatlab/        - toolbox by A. Banerjee and S. Sra used for sampling from von Mises-Fisher distribution
   |--- manopt/           - manifold optimization toolbox by N. Boumal and B. Mishra
   |--- matlab-xunit/     - xUnit tests for Matlab
   
 </pre>
## HOW TO CITE ##

If you are using this toolbox in your research please cite the following paper:
<pre>
@Article{mixest,
author = {Reshad Hosseini and Mohamadreza Mash'al},
title = {MixEst: An Estimation Toolbox for Mixture Models},
year = {to appear},
url = {http://visionlab.ut.ac.ir/mixest}
}
</pre>

Our toolbox depends on Manopt toolbox for doing manifold optimization, please cite the following paper if you are using manifold optimization methods of the toolbox.
<pre>
@Article{manopt,
  author  = {Nicolas Boumal and Bamdev Mishra and P.-A. Absil and Rodolphe Sepulchre},
  title   = {{M}anopt, a {M}atlab Toolbox for Optimization on Manifolds},
  journal = {Journal of Machine Learning Research},
  year    = {2014},
  volume  = {15},
  pages   = {1455--1459},
  url     = {http://www.manopt.org}
}
</pre>

If you are using manifold LBFGS method, please cite the following paper:
<pre>
@Article{gopt,
  title={Conic Geometric Optimization on the Manifold of Positive Definite Matrices},
  author={Sra, Suvrit and Hosseini, Reshad},
  journal={SIAM Journal on Optimization},
  volume={25},
  number={1},
  pages={713--739},
  year={2015},
  publisher={SIAM}
}
</pre>





## Overview ##

Using MixEst, you can fit arbitrary distributions (usually mixture models) to your data set by various optimization techniques, with minimal effort. The following is the list of what you can do with MixEst:

- Defining mixtures of arbitrary distributions
- Estimating the parameters of mixture models using EM algorithm or manifold optimization
- Avoiding local minima using split\&merge algorithms
- Simultaneous solving local minima and model selection using Competitive EM.
- Avoiding overfitting by penalized, stochastic manifold optimization or cross-validation
- Other tools like sampling, AIC, BIC, entropy, KL-divergence, etc.

<section>
<img src="http://visionlab.ut.ac.ir/mixest/docs/examples/img/example3.gif">
</section>
An example showing the toolbox fitting a mixture on von Mises-Fisher distributions.


<section>
<img src="http://visionlab.ut.ac.ir/mixest/docs/examples/img/example4.gif">
</section>
An example visualizing the estimation process of a modified version of the Competitive EM (CEM) algorithm on sample 2-D data.

-------------------------------------------------------------------------------
