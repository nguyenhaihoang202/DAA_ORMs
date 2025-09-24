# DAA_ORMs
This repository hosts the full codebase for my Master's thesis (in Bioinformatics and Computational Biology, University College Cork) on optimising Bayesian ordinal logistic regression models for microbiome differential-abundance analysis (DAA). The codebase implements two Bayesian ordinal models—BALOR (ordinal-only with Laplace shrinkage) and BZIOLR (zero-inflated extension)—and reproduces all simulation studies and real-data benchmarks (HMP 2012 gingival; Ravel 2011 BV).

This project adapted some code scripts from Juho Pelto (jepelt@utu.fi), and the adapted code scripts will be denoted with "JP".

The Stan models (BALOR and BZIOLR) are written with within-chain parallelisation enabled, so they can take advantage of multi-threading via CmdStanR. This allows faster sampling on modern multi-core CPUs, especially when analysing many taxa simultaneously. Thread usage can be controlled through the threads_per_chain argument in CmdStanR, and all simulation and benchmark scripts in this repository are pre-configured to exploit this option.

'Real_Benchmark' folder contains the codes used in benchmarking BALOR and BZIOLR with conventional DAA methods against datasets from MicrobiomeBenchmarkData (https://www.bioconductor.org/packages/release/data/experiment/html/MicrobiomeBenchmarkData.html)

The 'Simulation' folder contains all generators used in the thesis, the data are simulated from this study (https://doi.org/10.1038/nature24460), with 10 seeds to replicate. In this folder, the main scripts (denoted as 'Baseline') run five core methods side by side: frequentist ORM, ANCOM-BC, LinDA, BALOR, and BZIOLR. 

The 'Data Curation' and 'Model Development' folders are for some code scripts I worked on throughout the 4-month internship but were not included in the thesis.

The 'Figures' folder contains all the code scripts to generate the figures in the thesis.

Additional methods (e.g., DESeq2, MaAsLin2, corncob, LDM) are executed in separate scripts (denoted as 'Extra') so that they can be added without refitting the Bayesian models. A dedicated merge code ('Merge') then harmonises all results by normalising column formats (estimate, SE, p, q, CI bounds, significance) and combines baseline and extras into unified result objects.

The main results of the thesis can be replicated by using the following versions of the R packages: cmdstanr 0.9.0, posterior 1.6.1, tidyverse 2.0.0, rms 8.0-0 (for ORM), TreeSummarizedExperiment 2.14.0, phyloseq 1.50.0, ANCOMBC 2.8.1, LinDA 0.2.0, DESeq2 1.46.0, MaAsLin2 1.20.0, corncob 0.4.2, LDM 6.0.1. The codes ran in RStudio 2025.05.0, with R version 4.4.3.
