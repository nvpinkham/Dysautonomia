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


res <- NULL

for(i in 1 : length(treats)){
  for(j in 1 : i){
    if(i != j){

      par(mfrow=c(3,1))

      res.i <- pair.meta(map.meta, meta, comps = c(treats[i], treats[j]))

      res <- rbind(res, res.i)
    }
  }
}


res <- as.data.frame(res)
row.names(res) <- NULL
res <- res[order(res$group.1) , ]

write.csv(res, "Metabalome_between_groups.csv")
res <- read.csv("Metabalome_between_groups.csv")

res1 <- res[res$group.1 == "Control_Cohoused" & res$group.2 == "Mutant_Cohoused",]
res2 <- res[res$group.1 == "Control_Separate" & res$group.2 == "Mutant_Separate",]

res <- rbind(res1, res2)
res <- res[order(res$group.2),]

res$col <- as.numeric(as.factor(paste(res$group.1, res$group.2)))

colors <- colors()
colors <- colors[-grep("grey|gray", colors)]

set.seed(42)
res$col <- sample(colors)[res$col ]

res$n.scale <-   scale(res$n, center = T, scale = T)
res$n.scale <-   res$n.scale + diff(c(min(res$n.scale), 1))

res$dis <- as.factor(paste0(  res$group.1, "_",   res$group.2))
res$day <- res$day + (as.numeric(res$dis) * 10) - 15


#setEPS()       # Set postscript arguments
#postscript(paste0("PEMANOA_metabolite_summary",
#           width = 6, height = 8)


par(mfrow=c(2, 2))
par(mar = c(5.1, 4.1, 4.1, 2.1))


plot(res$p.val ~ res$day,
     pch = 21,
     cex = res$n.scale,
     bg = res$col,
     main = "Metabolite PERMANOVA",
     xlab = "Months ppost weaning",
     ylab = "P value",
     xlim= c(80, 360),
     ylim = c(-.1, .7),
     xaxt = "n")

axis(1, at = c(111, 201, 291),
     labels = c("3", "6", "9"))

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
                  mean(res$n.scale),
                  max(res$n.scale)),
       legend = paste("\n   n =", c(min(res$n),
                                    round(mean(res$n)),
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
     pch = 21,
     cex = res$n.scale,
     bg = res$col,
     main = "Metabolite PERMANOVA",
     xlab = "Months post weaning",
     ylab = "F stat",
     xlim = c(80, 360),
     ylim = c(.55, 2.95),
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
                  mean(res$n.scale),
                  max(res$n.scale)),
       legend = paste("\n   n =", c(min(res$n),
                                    round(mean(res$n)),
                                    max(res$n)),
                      "\n"),
       "right")

par(mar = c(0, 0, 0, 0))
plot.new()

res <- res[rev(order(res$p.val)), ]

# legend("left", bty="n", fill = unique(res$col),
# legend = unique(paste("\n", res$group.2, "\nto", res$group.1, "\n")))

table(map.meta$Genotype, map.meta$Age.Bin)













