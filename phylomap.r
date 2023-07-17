# plot phylogenetic tree on map
# reference: http://blog.phytools.org/2022/04/combining-contmap-and-phylotomap-plots.html

library(phytools)
library(ape)

setwd("...")

mytree <- read.tree("#.tre")
md <- read.csv("gis_data.csv")

rownames(md) <- md$X
head(md)
md1 <- md[,2:3]
head(md1)

obj <- phylo.to.map(mytree, md1, plot=FALSE)
obj

cols <- setNames(sample(rainbow(n=Ntip(mytree))),mytree$tip.label)
cols

plot(obj, colors = cols, ftype="i", fsize=0.6, cex.points=c(0.7,1.2))
