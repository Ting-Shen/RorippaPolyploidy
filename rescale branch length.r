library(ape)

setwd("...")
mytree1 <- read.tree("#.tre")

# drop outgroup

og <- c(...)

write.tree(drop.tip(mytree1,og),file = "#.tre", append = TRUE, digits = 10, tree.names = FALSE)

# rescaled trees for multiple times

# relaxed clock model

mytree <- read.tree("#.tre")
cal <- makeChronosCalib(mytree, node = "root", age.min = ..., age.max = ..., interactive = FALSE, soft.bounds = FALSE)

for (i in 1:1000) {
  rel <- chronos.control(nb.rate.cat = 10.0)
  chr.rel <- chronos(mytree, model = "relaxed", lambda = 1.0, calibration = cal, control = rel)
  a <- chr.rel
  myfile <- file.path("...", paste0("relaxed.", i,".tre"))
  write.tree(a, file = myfile)
}


# correlated
for (i in 1:1000) {
  chr.corr <- chronos(mytree, model = "correlated", lambda = 1.0, calibration = cal, control = chronos.control())
  a <- chr.corr
  myfile <- file.path("...", paste0("correlated.", i,".tre"))
  write.tree(a, file = myfile)
}

# discrete/ strict
for (i in 1:1000) {
  str <- chronos.control(nb.rate.cat = 1.0)
  chr.str <- chronos(mytree, model = "discrete", lambda = 1.0, calibration = cal, control = str)
  a <- chr.str
  myfile <- file.path("...", paste0("discrete.", i,".tre"))
  write.tree(a, file = myfile)
}

# creat and view an object with file names of *.

setwd("...")
f <- list.files(pattern = "*.tre")

# read all the list files and store them in one list named as d

d <- lapply(f, read.table)

# give names to the elements in the list

names(d) <- c(1:1000)

# export d

write.table(d, "rescaled_con.tre")

# violin plot of the loglik values to select the most fitted model

library(ggstatsplot)
library(palmerpenguins)
library(tidyverse)

md <- read.csv("rescale_result_loglik.csv", header = TRUE)

attach(md)
head(md)

plt <- ggbetweenstats(data = md, x=model, y=Loglik)
plt

plt <- plt +
  labs(
    x= "Model of substitution rate variation among branches",
    y= "Loglik"
  )

plt

