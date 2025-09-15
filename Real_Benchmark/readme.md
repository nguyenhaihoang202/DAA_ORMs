This folder contains the codes that work with the MicrobiomeBenchmarkData from GamboaTuz study.

Originally, the Bayesian models were saved under the names res_bayes4 (Laplace-ordinal) and res_zi (zero-inflated ordinal).

Later in the project, I decided to add additional methods (DESeq2, MaAsLin2, corncob, LDM) and to unify the naming convention across all scenarios. To avoid refitting the Bayesian models, I included small renaming helper blocks (rename_obj) inside the code. These blocks simply re-assign objects to the new names:

res_bayes4 → res_balor

res_zi → res_bziolr

This ensures consistency with the baseline/extra/merge workflow, where the five core methods are always labeled: ORM, ANCOM-BC, LinDA, BALOR, and BZIOLR.

By separating renaming from fitting, I can extend the analyses with extra methods and grid evaluations without re-running the expensive Bayesian sampling.
