# DAA_ORMs
This repository hosts the full codebase for my Master's thesis (in Bioinformatics and Computational Biology, University College Cork) on optimising Bayesian ordinal logistic regression models for microbiome differential-abundance analysis (DAA). The codebase implements two Bayesian ordinal models—BALOR (ordinal-only with Laplace shrinkage) and BZIOLR (zero-inflated extension)—and reproduces all simulation studies and real-data benchmarks (HMP 2012 gingival; Ravel 2011 BV).

Real_Benchmark folder contains includes the codes used in benchmarking BALOR and BZIOLR with conventional DAA methods against datasets from MicrobiomeBenchmarkData (https://www.bioconductor.org/packages/release/data/experiment/html/MicrobiomeBenchmarkData.html)

Simulations folder contains all generators used in the thesis: symmetric (no net load shift), unbalanced/realistic compositional effects, and single-taxon “bloomer,” plus niche scenarios for prevalence-only and partial global load-shift.

The main results of the thesis can be replicated by using the following versions of the R packages: cmdstanr 0.7.x (Stan 2.3x), posterior 1.5.0, tidyverse 2.0.0, rms 6.5-0 (for ORM), TreeSummarizedExperiment 2.10.0, phyloseq 1.44.0, ANCOMBC 2.8.1, LinDA 0.1.0, DESeq2 1.40.2, MaAsLin2 1.14.1, corncob 0.3.2, LDM 6.0.1.
