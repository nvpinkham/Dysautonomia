# Ex_fig_5
library(vegan)
library(labdsv)
library(vioplot)
library(cluster)
library(gplots)
library(RColorBrewer)
library(dendextend)
library(dplyr)

source("R/mice_functions.22.11.27.R")

'
Control Separate: #F0E442  (should be yellow color)
  FD Separate: #D55E00   (should be orange color)
  Control Cohoused: #CC79A7   (should be pink color)
  FD Cohoused: #009E73      (should be turquoise color)
'
# Mouse IDs 1600s and 1700s are first gen.
# Mouse IDs in the 4500s are second.
# Mouse IDs 0-300 are third.

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

for(i in 1:nrow(map.meta)){
  colfunc <- colorRampPalette(c("white",
                                map.meta$col[i],
                                "black"))
  col.pick <- colfunc(375)
  map.meta$Age.Col[i] <- col.pick[map.meta$Age.Bin[i]]

}

treats <-  c( "Mutant_Cohoused" , "Mutant_Separate",
              "Control_Cohoused", "Control_Separate")



map.1 <- get.close(map.meta, 100)
map.2 <- get.close(map.meta, 200)
map.3 <- get.close(map.meta, 300)

map.meta <- rbind(map.1,
                  map.2,
                  map.3)

res <-  nums <- NULL

map.meta$Age.Bin <- map.meta$Age.Bin - 21# DPW

par(mfrow=c(3,2))

aa <- colnames(meta)
bb <- substr(aa, 1, 1)
cc <- substr(aa, 2, 2)
bb[bb == "X" & grepl("[0-9]", cc)] <- ""
dd <- paste0(bb, substr(aa, 2, 90))
ee <- gsub("[.]", "-", dd)
colnames(meta) <- ee


colnames(meta)[3]

# setEPS()       # Set postscript arguments
# postscript(paste0("MC2MS_", Sys.Date(), ".eps"),width = 8, height = 10)
par(mfrow=c(3,2))

pair.meta3(map.meta, meta,
           comps = c("Mutant_Cohoused", "Mutant_Separate"))

# dev.off()


#setEPS()       # Set postscript arguments
#postscript(paste0("MC2CC_", Sys.Date(), ".eps"), width = 8, height = 10)

par(mfrow=c(3,2))
pair.meta3(map.meta, meta,
           comps = c("Mutant_Cohoused", "Control_Cohoused"))


# dev.off()

#setEPS()       # Set postscript arguments
#postscript(paste0("MS2CC_", Sys.Date(), ".eps"), width = 8, height = 10)

par(mfrow=c(3,2))
pair.meta3(map.meta, meta,
           comps = c("Control_Cohoused", "Mutant_Separate"))

# dev.off()




########## compile source data file ###########33

comps = c("Mutant_Cohoused", "Mutant_Separate")

dim(map.meta)
dim(meta)

meta <- meta[ row.names(meta) %in% map.meta$ID.Gen.Tmt.Age,]

all(map.meta$ID.Gen.Tmt.Age == row.names(meta))
pick <- map.meta$discription %in% comps

sup_data_8a <- cbind(map.meta$ID.Gen.Tmt.Age[pick],
                     map.meta$Genotype[pick],
                     map.meta$Treatment.type[pick],
                     map.meta$Age.Days[pick] - 21)

colnames(sup_data_8a) <- c("sample.id",
                           "genotype",
                           "housing",
                           "DPW")



sup_data_8a <- cbind("", "", sup_data_8a, as.matrix(meta)[pick,-1])
sup_data_8a[1,1] <- "Supplemental figure 8a"
sup_data_8a <- rbind(sup_data_8a, "", "")


write.csv(sup_data_8a,"sup_data_8a.csv", row.names = F)


######
########## compile source data file ###########33

comps = c("Mutant_Cohoused", "Control_Cohoused")
pick <- map.meta$discription %in% comps

sup_data_8b <- cbind(map.meta$ID.Gen.Tmt.Age[pick],
                     map.meta$Genotype[pick],
                     map.meta$Treatment.type[pick],
                     map.meta$Age.Days[pick] - 21)

colnames(sup_data_8b) <- c("sample.id",
                           "genotype",
                           "housing",
                           "DPW")

sup_data_8b <- cbind("", "", sup_data_8b, as.matrix(meta)[pick,-1])
sup_data_8b[1,1] <- "Supplemental figure 8b"

sup_data_8b <- rbind(sup_data_8b, "", "")

write.csv(sup_data_8b,"sup_data_8b.csv", row.names = F)

####
########## compile source data file ###########33

comps = c("Control_Cohoused", "Mutant_Separate")
pick <- map.meta$discription %in% comps

meta <- meta[ row.names(meta) %in% map.meta$ID.Gen.Tmt.Age,]

all(map.meta$ID.Gen.Tmt.Age == row.names(meta))

sup_data_8c <- cbind(map.meta$ID.Gen.Tmt.Age[pick],
                     map.meta$Genotype[pick],
                     map.meta$Treatment.type[pick],
                     map.meta$Age.Days[pick] - 21)

colnames(sup_data_8c) <- c("sample.id",
                           "genotype",
                           "housing",
                           "DPW")

sup_data_8c <- cbind("", "", sup_data_8c, as.matrix(meta)[pick,-1])
sup_data_8c[1,1] <- "Supplemental figure 8c"

sup_data_8c <- rbind(sup_data_8c, "", "")

write.csv(sup_data_8c,"sup_data_8c.csv", row.names = F)

system("cat sup_data_8* > sup_data_8.csv")
system("rm sup_data_8a.csv")
system("rm sup_data_8b.csv")
system("rm sup_data_8c.csv")



