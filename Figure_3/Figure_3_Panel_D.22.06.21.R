
library(vegan)
library(labdsv)
library(vioplot) # 0.3.7
library(cluster) # 2.1.3
library(gplots)
library(RColorBrewer)
library(dendextend)
library(dplyr)

source("R/mice_functions.22.06.22.R")

rs <-list.files("MetaboAnalystR-master/R", pattern = ".R", full.names = T)

for(i in rs[-1]){
 source(i)
}


dir.create("processing")

meta.6 <-  t(read.csv("data/6mo2018mousestool_umolepermg.csv",
                      header = T, row.names = 1))
meta.6 <- as.data.frame(meta.6)

meta.6 <- norm(row.names(meta.6), meta.6)

#########

p6 <- meta.6$Label

meta.6$Label <- NULL

adonis(dist(meta.6) ~ p6)

map.6 <- as.data.frame(cbind(row.names(meta.6),
                             p6,
                             as.numeric(as.factor(p6))))
colnames(map.6) <- c("id", "discription", "col")

map.6$col[map.6$discription == "FD_6mo"] <-  "#D55E00"
map.6$col[map.6$discription != "FD_6mo"] <-  "#F0E442"

plot.meta(meta.pick = meta.6, map.pick = map.6)

