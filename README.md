# Applying the Bayesian method DPGM-PN to the analysis of panel count data - Skintumor

## Overview
This repository contains R code implementing a novel Bayesian Markov Chain Monte Carlo (MCMC) method for the analysis of panel count data - Skintumor. Our approach is grounded in a Dirichlet Process Gamma Mixture Poisson Process model, aiming to provide robust insights into panel count data with an emphasis on accommodating within-subject correlation non-parametrically.

## Data
The primary dataset used to apply and validate our method is the skin tumor dataset, made available through the R package `PCDSpline`. This dataset originates from the paper by Yao, Bin, Lianming Wang, and Xin He, titled "Semiparametric regression analysis of panel count data allowing for within-subject correlation," published in Computational Statistics & Data Analysis, 97 (2016): 47-59.

## Methodology
Our Bayesian MCMC method leverages a Dirichlet Process Gamma Mixture Model to address a challenge inherent in panel count data analysis -- within-subject correlation. It generalizes the gamma frailty Poisson process model to allow an
unknown frailty distribution for analyzing panel count data. Specifically,
the frailty distribution is modeled nonparametrically by assigning a Dirichlet Process Gamma Mixture prior. An efficient Gibbs sampler is developed
to facilitate the Bayesian computation. This method stands out by allowing for a more flexible representation of the underlying distribution of the frailty, making it particularly suited for the complexities of panel count datasets.

### Key Features:
- **Dirichlet Process Gamma Mixture**: Enables a non-parametric approach to model frailty distributions, accommodating a wide range of variability and complexity within the data.
- **Poisson Process Model**: Ideal for handling panel count data, providing a natural framework for the analysis of events over time.
- **Bayesian MCMC Sampling**: The proposed Gibbs sampler works efficiently, which provides robust inference by quantifying uncertainty.

## Citation
Yao, Bin, Lianming Wang, and Xin He. "Semiparametric regression analysis of panel count data allowing for within-subject correlation." Computational Statistics & Data Analysis 97 (2016): 47-59.

## License
[MIT License](LICENSE.md)
