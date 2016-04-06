# STA-663-Final-Project
Spring 2016, Duke STA-663 Final Project
Gene Burin and Allen Ross
Paper: Account liability estimation improves power in ascertained case-control studies

## Abstract
The algorithm that will be implemented comes from the liability estimator as a phenotype (LEAP) framework that is utilized in genome-wide association studies (GWAS) which ensures an increase in power over the traditional linear mixed model (LMM) approach for case-control studies. This is of particular interest because GWAS have made a tremendous impact on understanding the role of genetics (SNPs) on phenotypes, but also present a multitude of opportunities to overcome statistical hurdles (e.g. type I error, power,  dependence/confounding, large dimensions). This algorithm addresses some of these issues, and builds off of LMMs, by using stochastic processes to estimate liability of a given phenotype. Although this is computationally intensive many parts are already performed by the regular LMM procedure, i.e. decomposition of the design matrix and the resulting eigendecomposition. Therefore, the power of the test of association can increase with relatively low added computational cost.

This project will provide a framework to implement parameter estimation, SVD, dimension reduction, and code optimization from an algorithmic standpoint. Accounting for inflation, confounding, and non-normality assumptions, particularly in rare diseases, to optimize statistical tests will be another focal point. Simulated data will be used to assess the improvement of this method over the standard LMM used in practice (power and speed). The development, implementation, and optimization of this algorithm will provide an important process in not only testing the statistical benefit, but also computational efficiency.    


## Background

## Implementation

## Optimization

## Results

## Comparisons

## References
[1] O. Weissbrod, C. Lippert, D. Geiger, D. Heckerman, *Accurate liability estimation improves power in ascertained case-control studies*, Nature Methods, 2015.  