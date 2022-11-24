
library(vegan)
library(labdsv)
library(vioplot) # 0.3.7
library(cluster) # 2.1.3
library(gplots)
library(RColorBrewer)
library(dendextend)
library(dplyr)

rs <-list.files("MetaboAnalystR-master/R", pattern = ".R", full.names = T)

source("R/mice_functions.22.09.20.R")

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

adonis2(dist(meta.6) ~ p6)

map.6 <- as.data.frame(cbind(row.names(meta.6),
                             p6,
                             as.numeric(as.factor(p6))))
colnames(map.6) <- c("id", "discription", "col")

map.6$col[map.6$discription == "FD_6mo"] <-  "#D55E00"
map.6$col[map.6$discription != "FD_6mo"] <-  "#F0E442"

plot.meta(meta.pick = meta.6, map.pick = map.6)

########## compile source data file ###########33
all(row.names(meta.6) == map.6$id)
map.6$genotype <- gsub("_6mo", "", map.6$discription)

source_data_3d <- cbind(map.6$id, map.6$genotype)
colnames(source_data_3d) <- c("sample.id", "genotype")
source_data_3d <- cbind("", "", source_data_3d, "", "")
source_data_3d[1,1] <- "figure 3d"

source_data_3d[1,6] <- "Metabolite concentrations"

source_data_3d <- cbind(source_data_3d, as.matrix(meta.6))

write.csv(source_data_3d,"source_data_3d.csv", row.names = F)

