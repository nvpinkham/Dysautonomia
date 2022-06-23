# Ex_fig_5

library(vegan)
library(labdsv)
library(vioplot)
library(cluster)
library(gplots)
library(RColorBrewer)
library(dendextend)
library(dplyr)

source("R/mice_functions.22.06.22.R")

rs <-list.files("MetaboAnalystR-master/R", pattern = ".R", full.names = T)

for(i in rs[-1]){
  source(i)
}



'
Control Separate: #F0E442  (should be yellow color)
  FD Separate: #D55E00   (should be orange color)
  Control Cohoused: #CC79A7   (should be pink color)
  FD Cohoused: #009E73      (should be turquoise color)
'
# Mouse IDs 1600s and 1700s are first gen.
# Mouse IDs in the 4500s are second.
# Mouse IDs 0-300 are third.
setwd("/Users/nickpinkham/Desktop/Dysautonomia/")

map <- read.csv("data/Map_CoHoused_mice.8.16.csv")
map$discription <- paste0(map$Genotype, "_", map$Treatment.type)
map$ID.Gen.Tmt.Age <- paste0("X", map$ID.Gen.Tmt.Age)
map$DOB <- as.Date(map$DOB , "%m/%d/%y")

map$col[map$Treatment.type =="Separate" & map$Genotype == "Mutant" ] <- '#D55E00'
map$col[map$Treatment.type =="Separate" & map$Genotype == "Control"] <- '#F0E442'
map$col[map$Treatment.type =="Cohoused" & map$Genotype == "Mutant"] <-  '#009E73'
map$col[map$Treatment.type =="Cohoused" & map$Genotype == "Control"] <- '#CC79A7'

metas <- list.files("data/2020 Mouse Metabolites norm by weight", pattern = "ALL", full.names = T)

a <- read.csv(metas[1], header = T, row.names = 1)
a <- as.data.frame(t(a))

b <- read.csv(metas[2], header = T, row.names = 1)
b <- as.data.frame(t(b))

c <- read.csv(metas[3], header = T, row.names = 1)
c <- as.data.frame(t(c))

colnames(b) <- colnames(a)
colnames(c) <- colnames(a)
meta <- as.data.frame(rbind(a, b, c))

rm(a)
rm(b)
rm(c)

map.meta <- map[match(row.names(meta), map$ID.Gen.Tmt.Age) , ]
all(map.meta$ID.Gen.Tmt.Age == row.names(meta))

map.meta$NMR.Date.Collected <- as.Date(map.meta$NMR.Date.Collected , "%m/%d/%y")
all(map.meta$ID.Gen.Tmt.Age == row.names(meta))


######## make unique col variation for each time point


treats <-  c( "Mutant_Cohoused" , "Mutant_Separate",
              "Control_Cohoused", "Control_Separate")


map.1 <- get.close(map.meta, 100)
map.2 <- get.close(map.meta, 200)
map.3 <- get.close(map.meta, 300)

map <- rbind(map.1,
             map.2,
             map.3)

res <-  nums <- NULL

map.meta$Age.Bin <- map.meta$Age.Bin - 21

treats

par(mfrow=c(3,2))

pair.meta3(map.meta, meta, 
           comps = c("Mutant_Cohoused", "Mutant_Separate"))


pair.meta3(map.meta, meta, 
           comps = c("Mutant_Cohoused", "Control_Cohoused"))


pair.meta3(map.meta, meta, 
           comps = c("Control_Cohoused", "Mutant_Separate"))

