# Software accompanying M-Y Chan's manuscript on "Improving EnsDA through PESE-GC"

PESE-GC stands for "Probit-space Ensemble Size Expansion".

EnsDA stands for "Ensemble Data Assimilation"

This manuscript will be submitted to Nonlinear Processes in Geophysics.

Author of software and manuscript: Man-Yau Chan (also known as "Joseph").


&nbsp; &nbsp;

## Manuscript abstract

Small forecast ensemble sizes (<100) are common in the ensemble data assimilation (EnsDA) component of geophysical forecast systems, thus limiting the error-constraining power of EnsDA. This study proposes an efficient and embarrassingly parallel method to generate additional ensemble members: the Probit-space Ensemble Size Expansion for Gaussian Copulas (PESE-GC; "peace gee see"). Such members are called "virtual members". PESE-GC utilizes the users' knowledge of the marginal distributions of forecast model variables. Virtual members can be generated from any (potentially non-Gaussian) multivariate forecast distribution that has a Gaussian copula. PESE-GC's impact on EnsDA is evaluated using the 40-variable Lorenz 1996 model, several EnsDA algorithms, several observation operators, a range of EnsDA cycling intervals and a range of forecast ensemble sizes. Significant improvements to EnsDA (p<0.01) are observed when either 1) the forecast ensemble size is small (<20 members), 2) the user selects marginal distributions that improves the forecast model variable statistics, and/or 3) the rank histogram filter is used with non-parametric priors in high forecast spread situations. These results motivate development and testing of PESE-GC for EnsDA with high-order geophysical models.



&nbsp; &nbsp;

## Directory structure
1) `manuscript_figs_2024` -- Contains Python scripts used to generate demonstrative and conceptual figures in the manuscript (Figs 2-6, and 10).
2) `lorenz96_expts` -- Contains a) scripts used to run the Lorenz 1996 experiments described in the manuscript, and b) scripts used to evaluate and visualize the performance of EnsDA with PESE-GC.
3) `DART_w_PESE-GC` -- NCAR Data Assimilation Research Testbed source code with PESE-GC implemented 


&nbsp; &nbsp;


## Instructions to replicate experiments in manuscript
See `lorenz96_expts/README.md`


&nbsp; &nbsp;

## Additional notes

PESE-GC is implemented in `DART_fortran_src/assimilation_code/modules/assimilation/pese_mod.f90`.

Some API-like documentation of PESE-GC is available in `DART_fortran_src/assimilation_code/modules/assimilation/pese_mod.rst`.
