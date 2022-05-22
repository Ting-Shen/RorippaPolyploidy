# Reconstructing Ancestral States for Discrete Variable (ploidy level, diploidy (state 1) or polyploidy (state 2)),
# under either the two ploidy levels model (TPL), and the multiple ploidy levels model (MPL).
# MPL model coupled with 3 hard-wired models: equal rate model (ER), all rates different model (ARD), and symmetric model (SYM).

# load libraries
library(ape)
library(phytools)

# set the work directory
setwd("...")

# read .tre file
mt <- read.tree("~.tre")

# read ploidy files
x2 <- read.csv("ploidy_2.csv", row.names = 1)
x5 <- read.csv("ploidy_5.csv", row.names = 1)
md2<-setNames(x2[,1],rownames(x2))
md5<-setNames(x5[,1],rownames(x5))

## for 2 ploidy levels
er2 <- ace(md2, mt, type = "discrete", model = "ER")
ard2 <- ace(md2, mt, type = "discrete", model = "ARD")
sym2 <- ace(md2, mt, type = "discrete", model = "SYM")

er2$loglik
ard2$loglik
sym2$loglik

# chi-square test
1-pchisq(2*abs(er2$loglik - ard2$loglik),1)
1-pchisq(2*abs(er2$loglik - sym2$loglik),1)
1-pchisq(2*abs(sym2$loglik - ard2$loglik),1)

# build a transition matrix for 2 ploidy levels, permitting the transition from diploidy (state 1) to polyploidy (state 2), but not the reverse
mytransition2 <- matrix(c(0,0,1,0), nrow = 2)
mytransition2

custome_asr2 <- ace(md2, mt, type = "discrete", model = mytransition2)
custome_asr2$loglik
custome_asr2$rates

# chi-square test
1-pchisq(2*abs(er2$loglik - custome_asr2$loglik),1)
1-pchisq(2*abs(ard2$loglik - custome_asr2$loglik),1)
1-pchisq(2*abs(sym2$loglik - custome_asr2$loglik),1)

# plot tree
plotTree(mt, fsize=0.6, ftype="i")

cols <- setNames(palette()[1:length(unique(md2))], sort(unique(md2)))
nodelabels(node = 1:mt$Nnode+Ntip(mt),
           pie = custome_asr2$lik.anc, piecol = cols, cex =0.5)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(mt)),fsize=0.8)

# simulating 1000 stochastic character
sim_asr2 <- make.simmap(mt, md2, model = mytransition2, nsim = 1000)
sim_asr2

pd_sim2 <- summary(sim_asr2, plot=FALSE)
pd_sim2
plot(pd_sim2, fsize=0.6, ftype="i")

# compare the posterior probabilities from stochastic mapping with our marginal ancestral states
plot(custome_asr2$lik.anc, pd_sim2$ace, xlab = "Marginal ancestral states",
     ylab = "Posterior probabilities from 100 stochastic mapping")

## for 5 ploidy levels
er5 <- ace(md5, mt, type = "discrete", model = "ER")
ard5 <- ace(md5, mt, type = "discrete", model = "ARD")
sym5 <- ace(md5, mt, type = "discrete", model = "SYM")

er5$loglik
ard5$loglik
sym5$loglik

# chi-square test
1-pchisq(2*abs(er5$loglik - ard5$loglik),1)
1-pchisq(2*abs(er5$loglik - sym5$loglik),1)
1-pchisq(2*abs(sym5$loglik - ard5$loglik),1)

# build a transition matrix for 5 ploidy levels, permitting the transition from low (state 1) to high (state 2), but not the reverse
# all equal rates
mytransition5.1 <- matrix(c(0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1,1,0), nrow = 5)
mytransition5.1

# all different rates
mytransition5.2 <- matrix(c(0,0,0,0,0,1,0,0,0,0,2,3,0,0,0,4,5,6,0,0,7,8,9,10,0), nrow = 5)
mytransition5.2

# equal rates from one ploidy level to others, but different rates for the other ploidy level to others
mytransition5.3 <- matrix(c(0,0,0,0,0,1,0,0,0,0,1,2,0,0,0,1,2,3,0,0,1,2,3,4,0), nrow = 5)
mytransition5.3

custome_asr5.1 <- ace(md5, mt, type = "discrete", model = mytransition5.1)
custome_asr5.1$loglik
custome_asr5.1$rates

custome_asr5.2 <- ace(md5, mt, type = "discrete", model = mytransition5.2)
custome_asr5.2$loglik
custome_asr5.2$rates

custome_asr5.3 <- ace(md5, mt, type = "discrete", model = mytransition5.3)
custome_asr5.3$loglik
custome_asr5.3$rates

# chi-square test
1-pchisq(2*abs(custome_asr5.1$loglik - custome_asr5.2$loglik),1)
1-pchisq(2*abs(custome_asr5.1$loglik - custome_asr5.3$loglik),1)
1-pchisq(2*abs(custome_asr5.2$loglik - custome_asr5.3$loglik),1)

1-pchisq(2*abs(er5$loglik - custome_asr5.1$loglik),1)
1-pchisq(2*abs(ard5$loglik - custome_asr5.1$loglik),1)
1-pchisq(2*abs(sym5$loglik - custome_asr5.1$loglik),1)

1-pchisq(2*abs(er5$loglik - custome_asr5.2$loglik),1)
1-pchisq(2*abs(ard5$loglik - custome_asr5.2$loglik),1)
1-pchisq(2*abs(sym5$loglik - custome_asr5.2$loglik),1)

1-pchisq(2*abs(er5$loglik - custome_asr5.3$loglik),1)
1-pchisq(2*abs(ard5$loglik - custome_asr5.3$loglik),1)
1-pchisq(2*abs(sym5$loglik - custome_asr5.3$loglik),1)

# plot tree
plotTree(mt, fsize=0.6, ftype="i")
cols <- setNames(palette()[1:length(unique(md5))], sort(unique(md5)))
nodelabels(node = 1:mt$Nnode+Ntip(mt),
           pie = custome_asr5.2$lik.anc, piecol = cols, cex =0.5)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(mt)),fsize=0.6)

plotTree(mt, fsize=0.6, ftype="i")
cols <- setNames(palette()[1:length(unique(md5))], sort(unique(md5)))
nodelabels(node = 1:mt$Nnode+Ntip(mt),
           pie = custome_asr5.1$lik.anc, piecol = cols, cex =0.5)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(mt)),fsize=0.6)

# simulating 1000 stochastic character under the most fitted model
sim_asr5.2 <- make.simmap(mt, md5, model = mytransition5.2, nsim = 1000)
sim_asr5.2

pd_sim5.2 <- summary(sim_asr5.2, plot=FALSE)
pd_sim5.2
plot(pd_sim5.2, fsize=0.6, ftype="i")

# compare the posterior probabilities from stochastic mapping with our marginal ancestral states
plot(custome_asr2$lik.anc, pd_sim2$ace, xlab = "Marginal ancestral states",
     ylab = "Posterior probabilities from 100 stochastic mapping")

# the end