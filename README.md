# STA-663-Final-Project
Spring 2016, Duke STA-663 Final Project
Gene Burin and Allen Ross
Paper: Accurate liability estimation improves power in ascertained case-control studies

## Abstract
The algorithm that will be implemented comes from the liability estimator as a phenotype (LEAP) framework that is utilized in genome-wide association studies (GWAS) which ensures an increase in power over the traditional linear mixed model (LMM) approach for case-control studies. This is of particular interest because GWAS have made a tremendous impact on understanding the role of genetics (SNPs) on phenotypes, but also present a multitude of opportunities to overcome statistical hurdles (e.g. type I error, power,  dependence/confounding, large dimensions). This algorithm addresses some of these issues, and builds off of LMMs, by using stochastic processes to estimate liability of a given phenotype. Although this is computationally intensive many parts are already performed by the regular LMM procedure, i.e. decomposition of the design matrix and the resulting eigendecomposition. Therefore, the power of the test of association can increase with relatively low added computational cost.

This project will provide a framework to implement parameter estimation, SVD, dimension reduction, and code optimization from an algorithmic standpoint. Accounting for inflation, confounding, and non-normality assumptions, particularly in rare diseases, to optimize statistical tests will be another focal point. Simulated data will be used to assess the improvement of this method over the standard LMM used in practice (power and speed). The development, implementation, and optimization of this algorithm will provide an important process in not only testing the statistical benefit, but also computational efficiency.    

## Background
Case-control genome-wide association studies (GWAS) focus on testing the association of single nucleotide polymorphisms (SNPs) and a given phenotype, which is binary. To increase the power of detecting an association it is beneficial to recode this phenotype as a continuous variable termed liability measure. The ability to accurately and efficiently estimate this liability measure is the purpose of this algorithm. There are many other frameworks to test association (linear regression, linear mixed model, etc.), but none are as robust, or lead to as much power, when dealing with confounding and other model assumptions as LEAP.   

The underlying estimation algorithm is only slightly more computationally intesive than any other methods listed above. LEAP is based off of a liability threshold model where the liability is estimated from a probit model, which is fitted using an iterative gradient descent algorithm with newton's method being used at each step. The genetic and environmental components of the liability measurement are generated from a joint MAP estimate, where only the genetic component is used in the algorithm to better adjust for related individuals. An estimate could be derived from the posterior mean for the liability measurement, but this is generally less accurate, and more computationally intensive. Runtime for original algorithm is approximately eight minutes per chromosome, insensitive to number of SNPs tested, for 8,000 individuals.  

## Implementation
Update code from Python 2.7 to be compatible with Python 3.xx


## References
[1] O. Weissbrod, C. Lippert, D. Geiger, D. Heckerman, *Accurate liability estimation improves power in ascertained case-control studies*, Nature Methods, 2015.  
[2] C. Lippert, J. Listgartern, Y. Liu, et al. *FaST linear mixed models for genome-wide association studies*, Nature Methods, 2011.
