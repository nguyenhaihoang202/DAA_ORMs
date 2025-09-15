# DAA_ORMs
This repository hosts the full codebase for my Master's thesis (in Bioinformatics and Computational Biology, University College Cork) on optimising Bayesian ordinal logistic regression models for microbiome differential-abundance analysis (DAA). The codebase implements two Bayesian ordinal models—BALOR (ordinal-only with Laplace shrinkage) and BZIOLR (zero-inflated extension)—and reproduces all simulation studies and real-data benchmarks (HMP 2012 gingival; Ravel 2011 BV).

The Stan models (BALOR and BZIOLR) are written with within-chain parallelisation enabled, so they can take advantage of multi-threading via CmdStanR. This allows faster sampling on modern multi-core CPUs, especially when analysing many taxa simultaneously. Thread usage can be controlled through the threads_per_chain argument in CmdStanR, and all simulation and benchmark scripts in this repository are pre-configured to exploit this option.

Real_Benchmark folder contains the codes used in benchmarking BALOR and BZIOLR with conventional DAA methods against datasets from MicrobiomeBenchmarkData (https://www.bioconductor.org/packages/release/data/experiment/html/MicrobiomeBenchmarkData.html)

Simulations folder contains all generators used in the thesis: symmetric (no net load shift), unbalanced/realistic compositional effects, and single-taxon “bloomer,” plus niche scenarios for prevalence-only and partial global load-shift.

The main results of the thesis can be replicated by using the following versions of the R packages: cmdstanr 0.9.0, posterior 1.6.1, tidyverse 2.0.0, rms 8.0-0 (for ORM), TreeSummarizedExperiment 2.14.0, phyloseq 1.50.0, ANCOMBC 2.8.1, LinDA 0.2.0, DESeq2 1.46.0, MaAsLin2 1.20.0, corncob 0.4.2, LDM 6.0.1. The codes ran in RStudio 2025.05.0, with R version 4.4.3.
