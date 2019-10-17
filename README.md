# Surface_approximation_GP

updated on October 17th, 2019

### The Matlab programs in this project

`Approx_Surface.m` is a draft program for any developing.

(Finished)

1. `Approx_Surface_GPR.m` is the very conventional Gaussian process regression program [Rasmussen 2003].

2. `Approx_Surface_poly55.m` uses the fitting toolbox of MATLAB for method comparison.

3. `Approx_Surface_GPR_datasetReduction.m` works with `datasetReduction.m`, which naively reduces the data by step.

4. `Approx_Surface_SVGP.m` is going to implement the SVGP method [Titsias 2009].

(On-going)

5. Modeling for experiment data.

### References:

Rasmussen C E. Gaussian processes in machine learning[C]//Summer School on Machine Learning. Springer, Berlin, Heidelberg, 2003: 63-71.

Burt D R, Rasmussen C E, Van Der Wilk M. Rates of Convergence for Sparse Variational Gaussian Process Regression[J]. arXiv preprint arXiv:1903.03571, 2019.

Titsias M. Variational learning of inducing variables in sparse Gaussian processes[C]//Artificial Intelligence and Statistics. 2009: 567-574.

### The Regression of Noisy Surface Data

<p float="left">
<img width="350" src=https://user-images.githubusercontent.com/26374671/61990819-6815b180-b082-11e9-93e7-20d8b26ffc7a.png alt="Noisy data" />
<img width="350" src=https://user-images.githubusercontent.com/26374671/61990831-95625f80-b082-11e9-93e5-d926ea108323.png alt="Original data" />
</p>
