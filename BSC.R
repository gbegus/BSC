bsc <- function (counts, surveyed, order = T, R = 10000) {
  # The function bsc() takes two vectors of equal length as arguments:
  # a vector with counts of languages with a sound changes required for
  # an alternation Ak, and a vector of languages surveyed for each sound change.
  # The function internally transforms the vectors with counts into a binomial 
  # distribution of successes and failures for each sound change in the count. 
  # It returns R  bootstrap replicates of the Historical Probability of A1,
  # computed according to Begus (2017). 
  # Stratified non-parametric bootstrapping is performed based on the boot package: 
  # the output of bsc() is an object of class 'boot'. The output of bsc() should be 
  # used as an argument of summary.bsc() (see summary.bsc()), which returns the observed
  # Px and 95% BC_a CIs. Two optional arguments of bsc() are order (if True, Historical
  # Probabilities are divided by n!) and R, which determines the number of bootstrap replicates.
  
  library(boot)
  if (length(counts) != length(surveyed)) {
    stop ("Vectors must be of equal length.")
  }
  binom <- unlist(mapply(c,
                         lapply(counts, function(x) rep(1, x)),
                         lapply(surveyed - counts, function(x) rep(0, x)), SIMPLIFY=F)
                  )
  
  snumb <- paste("s", 1:length(surveyed), sep = "")
  ident <- rep(snumb, surveyed)
  
  scsample <- data.frame(binom, ident)
  
  if (order == TRUE) {
    n <- factorial(length(counts))
  } else {
      n <- 1
      }
  
  bsc <- function(x, id) {
    sc1 <- tapply(x[id,1], x[id,2], mean)
    sc <- prod(sc1) / n
    return(sc)
  }
  
  boot.scsample <- boot(scsample, statistic = bsc, R, strata = scsample[, 2])
  return(boot.scsample)
}



bsc2 <- function(bsc.alt1a, bsc.alt2a, order = T, R = 10000){
  
  # The function bsc2() compares the Historical Probabilities 
  #of two processes with BSC. It takes as an input the output 
  #of bsc() for the process in question. The function transforms
  #the counts into a binomial distribution of successes and failures.
  #It returns R  bootstrap replicates of the difference in Historical
  #Probability between the two alternations, computed according to 
  #Begus (2017). Stratified non-parametric bootstrapping is performed 
  #based on the boot package: the output of bsc2() is an object of class `boot.
  #The output of bsc2() should be used as an argument of summary.bsc2()
  #(see summary.bsc2()), which returns the observed Px and 95% BCa CIs 
  #for the difference. If 95% BCa CIs fall above or below zero, 
  #it spells out that the difference is significant, 
  #and that it is not otherwise.  Two optional arguments of bsc()
  #are order (if True, Historical Probabilities are divided by n!)
  #and R, which determines the number of bootstrap replicates.
  
  library(boot)
  bsc.alt1 <- bsc.alt1a$data
  bsc.alt2 <- bsc.alt2a$data
  bsc.alt1$scid <- "first"
  bsc.alt2$scid <- "second"
  bsc.diff.df <- rbind(bsc.alt1, bsc.alt2)
  bsc.diff.df$comb <- as.factor(paste(bsc.diff.df$scid, bsc.diff.df$ident, sep = ""))
  
  bsc.diff.df$scid <- NULL
  bsc.diff.df$ident <- NULL
  
  if (order == TRUE) {
    n1 <- factorial(length(unique(bsc.alt1$ident)))
    n2 <- factorial(length(unique(bsc.alt2$ident)))
    } else {
    n1 <- 1
    n2 <- 1
    }
  
  l <- length(unique(bsc.alt1$ident))
  m <- length(unique(bsc.alt2$ident))

  bsc.diff <- function(x, id) {
    sc1 <- tapply(x[id,1], x[id,2], mean)
    sca <- (prod(sc1[1:l]) / n1)
    scb <- (prod(sc1[(l+1):(l+m)]) / n2)
    sc <- sca - scb
    return(sc)
  }

  boot.diff <- boot(bsc.diff.df, statistic = bsc.diff, R, strata = bsc.diff.df[, 2])
  return(boot.diff)
}


summary.bsc <- function (bsc.alt) {
  
  # The function summary.bsc() computes the 95% BCa CI for the bootstrap
  # replicates based on the bsc() function (see bsc()) using the boot.ci()
  # function from the boot package and returns the observed and estimated 
  # Historical Probabilities.
  
  bsc.ci.alt <- boot.ci(bsc.alt, type="bca")
  title <- "BOOTSTRAPPING SOUND CHANGES"
  prob <- paste("Observed P =", round(bsc.alt$t0*100, digits = 5), "%")
  bca <- paste("Estimated 95 % BCa CI = [", round(bsc.ci.alt$bca[4]*100, digits = 4),"%,",
               round(bsc.ci.alt$bca[5]*100, digits = 4),"%]")
  probbca <- paste(prob, bca, sep = "\n")
  cat(title, probbca, sep = "\n\n")
}


summary.bsc2 <- function (bsc2.alt) {
  
  # The function summary.bsc2() computes the 95% BCa CI for the bootstrap
  # replicates based on the bsc2() function (see bsc2()) using the boot.ci()
  # function from the boot package and returns the observed and estimated differences
  # in Historical Probabilities of two alternations. 
  
  bsc2.ci.alt <- boot.ci(bsc2.alt, type="bca")
  title <- "BOOTSTRAPPING SOUND CHANGES - COMPARE"
  prob <- paste("Observed Delta P =", round(bsc2.alt$t0*100, digits = 5), "%")
  bca <- paste("Estimated 95 % BCa CI = [", round(bsc2.ci.alt$bca[4]*100, digits = 4),"%,",
               round(bsc2.ci.alt$bca[5]*100, digits = 4),"%]")
  
  if (bsc2.ci.alt$bca[4] > 0 & bsc2.ci.alt$bca[5] > 0) {
    sig <- "P(A1) is significantly higher than P(A2)."
  }
  else if (bsc2.ci.alt$bca[4] < 0 & bsc2.ci.alt$bca[5] < 0) {
    sig <- "P(A1) is significantly lower than P(A2)."
  } else {
    sig <- "P(A1) and P(A2) are not significantly different."
  }
  
  probbca <- paste(prob, bca, sep = "\n")
  cat(title, probbca,sig, sep = "\n\n")
}


plot.bsc <- function (bsc.alt, Alternation = c("Alternation")) {
  
  # The function plot.bsc() takes the output of bsc() as  input and 
  # plots the distribution of bootstrap replicates with the observed
  # Historical Probability of the process (solid line) and 95% BCa CI (dashed line),
  # calculated with the boot.ci() function from the boot package. 
  # The plotting is based on the ggplot2 package (Wickham 2009). 
  # An optional argument Alternation allows for the change of 
  # the name of the alternation in the legend.
  
  library(ggplot2)
  bsc.ci.alt <- boot.ci(bsc.alt, type = "bca")
  bsc.alt.df <- data.frame(bsc.alt$t)
  bsc.alt.df$name <- Alternation
  names(bsc.alt.df) <- c("boot", "Alternation")
  boot.plot <- ggplot(bsc.alt.df, aes(boot, fill = Alternation)) + geom_density(alpha = 0.5)  +
    geom_vline(xintercept = bsc.alt$t0, colour="red",linetype = "solid") +
    geom_vline(xintercept = bsc.ci.alt$bca[4], 
               colour="red", linetype = "dashed") +
    geom_vline(xintercept = bsc.ci.alt$bca[5], 
               colour="red", linetype = "dashed") +
    theme_bw() + xlab("Px in %") + ylab("") 
  return(boot.plot)
}



plot.bsc2 <- function (bsc.alt1, bsc.alt2, Alternation = c("Alternation 1","Alternation 2")) {
  
  # The function plot.bsc2() takes the output of bsc() as its input
  # (two alternations) and plots the distribution of bootstrap replicates 
  # with the observed Historical Probability of the process (solid line) 
  # and 95% BCa CI (dashed line), calculated with the boot.ci() function
  # from the boot package for each alternation. The plotting is based on 
  # the ggplot2 package (Wickham 2009). An optional argument Alternation
  # allows for the change of the name of the two alternations in the legend. 
  # Note that plot.bsc2() does not plot bootstrap replicates of the difference
  # between two Historical Probabilities, but rather bootstrap replicates of 
  # Historical Probabilities of each of the two alternations. To plot the bootstrap
  # replicates of the difference between two Historical Probabilities, apply plot.bsc()
  # to the output  of bsc2().
  
  library(ggplot2)
  bsc.ci.alt1 <- boot.ci(bsc.alt1, type = "bca")
  bsc.ci.alt2 <- boot.ci(bsc.alt2, type = "bca")
  bsc.alt1.df <- data.frame(bsc.alt1$t)
  bsc.alt2.df <- data.frame(bsc.alt2$t)
  bsc.alt1.df$name <- Alternation[1]
  bsc.alt2.df$name <- Alternation[2]
  names(bsc.alt1.df) <- c("boot", "Alternation")
  names(bsc.alt2.df) <- c("boot", "Alternation")
  bsc.alt.df<-rbind(bsc.alt1.df, bsc.alt2.df)
  boot.plot <- ggplot(bsc.alt.df, aes(boot, fill = Alternation)) +
    geom_density(alpha = 0.5)  +
    geom_vline(xintercept = bsc.ci.alt1$bca[4], 
               colour = "red", linetype = "dashed") +
    geom_vline(xintercept = bsc.ci.alt1$bca[5], 
               colour = "red", linetype = "dashed") +
    geom_vline(xintercept = bsc.alt1$t0, 
               colour = "red", linetype = "solid") +
    geom_vline(xintercept = bsc.ci.alt2$bca[4], 
               colour = "turquoise4", linetype = "dashed") +
    geom_vline(xintercept = bsc.ci.alt2$bca[5], 
               colour = "turquoise4", linetype = "dashed") +
    geom_vline(xintercept = bsc.alt2$t0, 
               colour = "turquoise4", linetype = "solid") +
    theme_bw() + xlab("Px in %") + ylab("") 
  
  return(boot.plot)
}



# Example:
pnd.counts <- c(97, 18, 27)
pnd.surveyed <- c(294, 263, 216)

pnd <- bsc(pnd.counts, pnd.surveyed)
summary.bsc(pnd)

# Output:
##BOOTSTRAPPING SOUND CHANGES
##
##Observed P = 0.04704 %
##Estimated 95 % BCa CI = [ 0.0261 %, 0.0862 %]


# Example:
pnd.counts <- c(97,18,27)
pnd.surveyed <- c(294,263,216)

fv.counts <- c(6,32,27)
fv.surveyed <- c(294,294,88)

pnd <- bsc(pnd.counts, pnd.surveyed)
fv <- bsc(fv.counts, fv.surveyed)

pndfv <- bsc2(pnd, fv)
summary.bsc2(pndfv)

# Output:
##BOOTSTRAPPING SOUND CHANGES - COMPARE
##
##Observed Delta P = 0.03568 %
##Estimated 95 % BCa CI = [ 0.0114 %, 0.0744 %]
##
##P(A1) is significantly higher than P(A2).


# Example:
pnd.counts <- c(97, 18, 27)
pnd.surveyed <- c(294, 263, 216)

pnd <- bsc(pnd.counts, pnd.surveyed)
plot.bsc(pnd, alternation = "PND")


#Example:
pnd.counts <- c(97,18,27)
pnd.surveyed <-c (294,263,216)

fv.counts <- c(6,32,27)
fv.surveyed <- c(294,294,88)

pnd <- bsc(pnd.counts, pnd.surveyed)
fv <- bsc(fv.counts, fv.surveyed)

plot.bsc2(pndfv)