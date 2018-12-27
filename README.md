# Replication Code: quantile_coherency_replication_pack

R code to replicate the Figures from the Barunik and Kley (2018):
*Quantile Coherency: A General Measure for Dependence between Cyclical
Economic Variables*, forthcoming in The Econometrics Journal, DOI:

The manuscript
[available here for download](https://ideas.repec.org/p/arx/papers/1510.06946.html)

## Requirements

To replicate the results, a number of other R packages have to be installed.
The script *install_required_packages.R* can be used to prepare for running
the replication scripts.

For replication of the Figure 4 showing cross-quantilogram, a package that 
is not available via the [CRAN](https://cran.r-project.org/mirrors.html).
Instead it is available from the authors homepage and needs to be
[downloaded](https://sites.google.com/site/whangyjhomepage/research/software)
and installed prior to running *install_required_packages.R*.

Replication of Figure 3 requires 32GB of memory to be available. All other
figures can be replicated with 16GB of memory being available. 

## Content

Note that the files
*Code_to_Figure_1.R*,
*Code_to_Figure_2.R*,
*Code_to_Figure_3.R*,
*Code_to_Figure_4.R*,
*Code_to_Figure_5.R*, and
*Code_to_Figures_6to8.R* can be used to replicate
the Figures 1, 2, 3, 4, 5 and Figures 6-8 respectively,
while files
*Code_to_Figure_S1.R*,
*Code_to_Figure_S2.R*,
*Code_to_Figure_S3.R*, and
*Code_to_Figure_S4.R* replicate
the Figures 1, 2, 3, and 4 in the Online Supplement.
