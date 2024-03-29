# BSC

Code for Beguš, Gašper. 2020. [Estimating historical probabilities of natural and unnatural processes](https://www.cambridge.org/core/journals/phonology/article/estimating-historical-probabilities-of-natural-and-unnatural-processes/CAC83CE585B82836CBC0D3A4BCDF17EB#article). *Phonology* 37(4): 515-549. doi:10.1017/S0952675720000263 

### bsc()

The function `bsc()` takes two vectors of equal length as arguments: a vector with counts of languages with a sound changes required for an alternation A<sub>k</sub>, and a vector of languages surveyed for each sound change. The function internally transforms the vectors with counts into a binomial distribution of successes and failures for each sound change in the count. It returns R  bootstrap replicates of the Historical Probability of A<sub>1</sub>, computed according to Begus (2019). Stratified non-parametric bootstrapping is performed based on the boot package: the output of `bsc()` is an object of class 'boot'. The output of `bsc()` should be used as an argument of `summary.bsc()`, which returns the observed P<sub>x</sub> and 95% BC<sub>a</sub> CIs. Two optional arguments of `bsc()` are order (if True, Historical Probabilities are divided by n!) and R, which determines the number of bootstrap replicates.


The function summary.bsc() computes the 95% BCa CI for the bootstrap replicates based on the `bsc()` function (see `bsc()`) using the `boot.ci()` function from the boot package and returns the observed and estimated Historical Probabilities.

Example: 
```
pnd.counts <- c(47,15,17)
pnd.surveyed <-c (294,263,216)

pnd <- bsc(pnd.counts, pnd.surveyed)
summary.bsc(pnd)

> BOOTSTRAPPING SOUND CHANGES
>
> Observed P = 0.01196 %
> Estimated 95 % BCa CI = [ 0.0059 %, 0.025 %]
```

### bsc2()

The function `bsc2()` compares the Historical Probabilities of two processes with BSC. It takes as an input the output of `bsc()` for the process in question. The function transforms the counts into a binomial distribution of successes and failures. It returns R bootstrap replicates of the difference in Historical Probability between the two alternations, computed according to Begus (2019). Stratified non-parametric bootstrapping is performed based on the boot package: the output of `bsc2()` is an object of class 'boot'. The output of `bsc2()` should be used as an argument of `summary.bsc2()`, which returns the observed Px and 95% BC<sub>a</sub> CIs for the difference. If 95% BCa CIs fall above or below zero, it spells out that the difference is significant, and that it is not otherwise. Two optional arguments of `bsc()` are order (if True, Historical Probabilities are divided by n!) and R, which determines the number of bootstrap replicates.

 The function `summary.bsc2()` computes the 95% BC<sub>a</sub>  CI for the bootstrap replicates based on the `bsc2()` function using the `boot.ci()` function from the *boot* package and returns the observed and estimated differences in Historical Probabilities of two alternations.

Example: 
```
pnv.counts <- c(28)
pnv.surveyed <- c(294)

pnv <- bsc(pnv.counts, pnv.surveyed)
summary.bsc(pnv)

pnd.counts <- c(47,15,17)
pnd.surveyed <-c (294,263,216)

pnd <- bsc(pnd.counts, pnd.surveyed)
summary.bsc(pnd)

pnvpnd <- bsc2(pnv, pnd)
summary.bsc2(pnvpnd)

> BOOTSTRAPPING SOUND CHANGES - COMPARE
> 
> Observed Delta P = 9.51185 %
> Estimated 95 % BCa CI = [ 6.4504 %, 13.2508 %]
>
> P(A1) is significantly higher than P(A2).

```

### plot.bsc()
The function `plot.bsc()` takes the output of `bsc()` as  input and plots the distribution of bootstrap replicates with the observed Historical Probability of the process (solid line) and 95% BC<sub>a</sub>  CI (dashed line), calculated with the boot.ci() function from the boot package. The plotting is based on the *ggplot2* package (Wickham 2009).  An optional argument Alternation allows for the change of  the name of the alternation in the legend.

```
plot.bsc(pnd)
```
![Plot](/pndBscPlot.jpeg)


### plot.bsc2()

 The function `plot.bsc2()` takes the output of `bsc()` as its input (two alternations) and plots the distribution of bootstrap replicates with the observed Historical Probability of the process (solid line) and 95% BC<sub>a</sub> CI (dashed line), calculated with the `boot.ci()` function from the *boot* package for each alternation. The plotting is based on the *ggplot2* package (Wickham 2009). An optional argument Alternation allows for the change of the name of the two alternations in the legend. Note that `plot.bsc2()` does not plot bootstrap replicates of the difference between two Historical Probabilities, but rather bootstrap replicates of Historical Probabilities of each of the two alternations. To plot the bootstrap replicates of the difference between two Historical Probabilities, apply `plot.bsc()` to the output  of bsc2().
 
 
 
```
ivv.counts <- c(38)
ivv.surveyed <- c(294)

ivv <- bsc(ivv.counts, ivv.surveyed)
summary.bsc(ivv)

plot.bsc2(pnv,ivv)
```
 
 ![Plot2](/pndBscPlot2.jpeg)
 

