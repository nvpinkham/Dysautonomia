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

map$X16SFile.ID[map$X16SFile.ID == "" ] <- NA

map$col[map$Treatment.type =="Separate" & map$Genotype == "Mutant" ] <- '#D55E00'
map$col[map$Treatment.type =="Separate" & map$Genotype == "Control"] <- '#F0E442'
map$col[map$Treatment.type =="Cohoused" & map$Genotype == "Mutant"] <-  '#009E73'
map$col[map$Treatment.type =="Cohoused" & map$Genotype == "Control"] <- '#CC79A7'

otu <- read.csv("data/FD_OTU_mice21.csv",
                header = T, row.names = 1)

not.in.map <- otu[!row.names(otu) %in% map$ID.Gen.Tmt.Age,]
row.names(not.in.map) # these mice are not from the third generation
rm(not.in.map)

map.otu <- map[!is.na(map$X16SFile.ID) , ]
map.otu <- map.otu[map.otu$DOB > "2018-12-31" , ] # remove 2nd gen mice!
#  Mouse IDs 0-300 are third.

otu <- otu[match(map.otu$ID.Gen.Tmt.Age, row.names(otu)) , ]
map.otu <- map.otu[match(map.otu$ID.Gen.Tmt.Age, row.names(otu)) , ]

all(map.otu$ID.Gen.Tmt.Age ==  row.names(otu))

map.otu$invsimp <- diversity(otu, "invsimp")

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
for(i in 1:nrow(map.otu)){
        colfunc <- colorRampPalette(c("white",
                                      map.otu$col[i],
                                      "black"))
        col.pick <- colfunc(375)
        map.otu$Age.Col[i] <- col.pick[map.otu$Age.Bin[i]]
}

for(i in 1:nrow(map.meta)){
        colfunc <- colorRampPalette(c("white",
                                      map.meta$col[i],
                                      "black"))
        col.pick <- colfunc(375)
        map.meta$Age.Col[i] <- col.pick[map.meta$Age.Bin[i]]

}

treats <-  c( "Mutant_Cohoused" , "Mutant_Separate",
              "Control_Cohoused", "Control_Separate")


map.1 <- get.close(map.otu, 100)
map.2 <- get.close(map.otu, 200)
map.3 <- get.close(map.otu, 300)



map <- rbind(map.1,
             map.2,
             map.3)

res <- NULL

for(i in 1 : length(treats)){
        for(j in 1 : i){
                if(i != j){

                        res.i <- pair.otu2(map, otu,
                                           c(treats[i], treats[j]))
                        res <- rbind(res, res.i)
                }
        }
}

res <- as.data.frame(res)
row.names(res) <- NULL
res <- res[order(res$DPW) , ]
write.csv(res, "Microbiome_between_groups.csv")

#setEPS()       # Set postscript arguments
#postscript(paste0("PERNOVA_16S_summary",
#           width = 6, height = 8)

res <- read.csv("Microbiome_between_groups.csv")

colnames(res)[2:3] <- c("group.1", "group.2")

res1 <- res[res$group.1 == "Control_Cohoused" & res$group.2 == "Mutant_Cohoused",]

res2 <- res[res$group.1 == "Control_Separate" & res$group.2 == "Mutant_Separate",]


res <- rbind(res1, res2)
res <- res[order(res$n.group.2),]

res$col <- as.numeric(as.factor(paste(res$group.1, res$group.2)))

colors <- colors()
colors <- colors[-grep("grey|gray", colors)]

set.seed(42)
res$col <- sample(colors)[res$col ]

res$n.scale <-   scale(res$n, center = T, scale = T)
res$n.scale <-   res$n.scale + diff(c(min(res$n.scale), 1))

res$dis <- as.factor(paste0(  res$group.1, "_",   res$group.2))

res$day <- res$DPW + (as.numeric(res$dis) * 10) - 15

par(mfrow=c(2, 2))
par(mar = c(5.1, 4.1, 4.1, 2.1))

#plot(res$p.val ~ res$day,
plot(res$p.val ~ res$day,
     xlim = c(80, 380),
     ylim = c(-.01, .1),
     pch = 21,
     cex = res$n.scale,
     bg = res$col,
     main = "16S PERMANOVA",
     xlab = "DPW",
     ylab = "P value",
     xaxt = "n")

axis(1, at = c(100, 200, 300),
     labels = c("79", "179", "279"))

comps <- paste(res$group.1, res$group.2)

abline(h = 0.05, col = 2, lty = 2)

for(i in 1 : 2){

        res.i <- res[comps == unique(comps)[i] , ]
        points(  res.i $p.val ~ res.i$day, bg = res.i$col[1], pch = 21,  cex = res.i$n.scale)
        points(  res.i $p.val ~ res.i$day, col = res.i$col[1], type = "l")
}

legend(pch = 21, pt.bg  = 8,
       bty="n", cex = .75,
       pt.cex = c(min(res$n.scale),
                  median(res$n.scale),
                  max(res$n.scale)),
       legend = paste("\n   n =", c(min(res$n),
                                    median(res$n),
                                    max(res$n)),
                      "\n"),
       "topright")

par(mar = c(0, 0, 0, 0))
plot.new()

#res <- res[rev(order(res$p.val)), ]

legend("left", bty="n",
       pch = 21, pt.bg  = unique(res$col),
       legend = c("FD vs. control mice cohoused",
                  "FD vs. control mice housed seperately"))
#######################################################

par(mar = c(5.1, 4.1, 4.1, 2.1))

plot(res$f.stat ~ res$day,
     xlim = c(80, 380),
     ylim = c(1.5, 5.5),
     pch = 21,
     cex = res$n.scale,
     bg = res$col,
     main = "16S PERMANOVA",
     xlab = "DPW",
     ylab = "F stat",
     xaxt = "n")


axis(1, at = c(100, 200, 300),
     labels = c("79", "179", "279"))

#axis(1, at = c(111, 201, 291),
#     labels = c("3", "6", "9"))

for(i in 1 : 2){

        res.i <- res[comps == unique(comps)[i] , ]
        points(  res.i $f.stat ~ res.i$day, bg = res.i$col[1], pch = 21,  cex = res.i$n.scale)
        points(  res.i $f.stat ~ res.i$day, col = res.i$col[1], type = "l")
}

legend(pch = 21, pt.bg  = 8,
       bty="n", cex = .75,
       pt.cex = c(min(res$n.scale),
                  median(res$n.scale),
                  max(res$n.scale)),
       legend = paste("\n   n =", c(min(res$n),
                                    median(res$n),
                                    max(res$n)),
                      "\n"),
       "right")


table(map.otu$Treatment.type, map.otu$Genotype, map.otu$Age.Bin)
table(map.otu$Genotype, map.otu$Age.Bin)





