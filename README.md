# BSC

Code for Begus, Gasper. 2019. [*Estimating historical probabilities of natural and unnatural processes*](https://ling.auf.net/lingbuzz/004299). Ms, UC Berkeley.

The function `bsc()` takes two vectors of equal length as arguments: a vector with counts of languages with a sound changes required for an alternation A<sub>k</sub>, and a vector of languages surveyed for each sound change. The function internally transforms the vectors with counts into a binomial distribution of successes and failures for each sound change in the count. It returns R  bootstrap replicates of the Historical Probability of A<sub>1</sub>, computed according to Begus (2019). Stratified non-parametric bootstrapping is performed based on the boot package: the output of `bsc()` is an object of class 'boot'. The output of `bsc()` should be used as an argument of `summary.bsc()` (see `summary.bsc()`), which returns the observed P_x and 95% BC<sub>a</sub> CIs. Two optional arguments of `bsc()` are order (if True, Historical Probabilities are divided by n!) and R, which determines the number of bootstrap replicates.
