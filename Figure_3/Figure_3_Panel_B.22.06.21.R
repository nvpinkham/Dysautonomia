
library(vegan)

map <- read.csv("data/ASD_map.csv",
                header=TRUE,
                row.names=1)

map$sample.id <- paste0(map$mouse,
                        ".",
                        substr(map$genotype, 1, 1),
                        ".",
                        map$time)

otu <- read.csv("data/FD_OTU_mice21.csv",
                row.names = 1,
                header = T)


otu <- otu[match(map$sample.id, row.names(otu)),]


map$cage <- paste0(map$genotype, ".", map$sex)

map$invsimp <- diversity(otu, "invsimp")
m.dist <- vegdist(otu)
map$discription <- paste("day", map$time, map$genotype)

map$color[map$genotype == "mutant"] <-  "#D55E00"
map$color[map$genotype != "mutant"] <-  "#F0E442"

mice <- unique(map$mouse)

map$bc.from.start <- NA
map$intial.invsimpt <- NA


for(i in 1 : length(mice)){

  otu.pick <- otu[map$mouse == mice[i],]
  map.pick <- map[map$mouse == mice[i],]

  res <- as.matrix(vegdist(otu.pick))[ 1 , ]

  aa <- map.pick$invsimp[map.pick$time == 0]
  if(length(aa)== 1){

    for(j in 1 : length( res)){
      # add distance tp map one at a time
      map$bc.from.start[row.names(map) == names(res)[j]] <-  res [j]

    }

    map$intial.invsimp[map$mouse == mice[i]] <-
      map.pick$invsimp[map.pick$time == 0]
  }
}







asd.summary <- matrix(nrow = 5, ncol = 15)
colnames(asd.summary) <- c("PERMANOVA p val", "F stat", "PERMANOVA DF",
                           "bc mu to WT", "bc mu",
                           "bc WT", "mean bc",
                           "mu.centroid", "control.centroid", "bc.centroids",
                           "wilk.mu", "wilk.control",
                           "t test p val", "t stat", "t DF")
rownames(asd.summary) <- unique(map$time)
asd.summary <-  as.data.frame(asd.summary)

time.line <- sort(unique(map$time))
dist.list <- list()
bb <- list()

for(i in 1 : 5){
  # summarize each timepoint
  now <- time.line[i]

  otu.pick <- otu[map$time == now , ]
  map.pick <- map[map$time == now , ]

  assign(paste0("otu.", now), otu.pick)
  assign(paste0("map.", now), map.pick)

  # PERMANOVA
  res <- adonis2(otu.pick ~ map.pick$genotype, permutations = 9999)
  asd.summary $`PERMANOVA p val`[i] <- res$`Pr(>F)`[1]
  asd.summary $`F stat`[i] <- res$F[1]
  asd.summary $`PERMANOVA DF`[i] <- res$Df[3]

  # explore distance matrix
  pick.dist <- as.matrix(vegdist(otu.pick))
  between.dist <- pick.dist[map.pick$genotype == "control",
                            map.pick$genotype ==  "mutant"]

  mu.dist <- as.dist(pick.dist[map.pick$genotype == "mutant",
                               map.pick$genotype ==  "mutant"])

  con.dist <- as.dist(pick.dist[map.pick$genotype == "control",
                                map.pick$genotype ==  "control"])


  asd.summary$`bc mu to WT`[i] <- mean(between.dist)
  asd.summary$`bc mu`[i] <- mean(mu.dist)
  asd.summary$`bc WT`[i] <- mean(con.dist)
  asd.summary$bc[i] <- mean(pick.dist)# overall BC

  # evaluate beta dispersion
  beta.res <- betadisper(vegdist(otu.pick), map.pick$genotype, type = "centroid")

  if(all(names(beta.res$distances) == row.names(map.pick))){
    map$bc.centroid[row.names(map) %in% names(beta.res$distances)] <- beta.res$distances
  }

  asd.summary$mu.centroid[i] <- mean(beta.res$distances[beta.res$group == "mutant"])
  asd.summary$control.centroid[i] <- mean(beta.res$distances[beta.res$group == "control"])
  asd.summary$bc.centroids[i] <- dist(beta.res$centroids)


  # perform T test; is their a difference in how deterministic/selective
  # the are the mice are of their microbiomes.

  # qqnorm(con.dist)# all looked normal
  # qqnorm(mu.dist)

  asd.summary$wilk.mu[i] <- shapiro.test(mu.dist)$p.value# used to reject NULL hypothesis of normality
  asd.summary$wilk.control[i] <- shapiro.test(con.dist)$p.value



  t.res <- t.test(as.vector(mu.dist),  as.vector(con.dist))

  asd.summary$`t test p val`[i] <- t.res$p.value
  asd.summary$`t stat`[i] <- t.res$statistic
  asd.summary$`t DF`[i] <- t.res$parameter

  w.res <- wilcox.test(as.vector(mu.dist),  as.vector(con.dist))
  asd.summary$`w test p val`[i] <- w.res$p.value

  # Save lists of distance values for box/violin plots

  aa  <- list(as.vector(mu.dist),
              as.vector(con.dist),
              as.vector(between.dist))

  names(aa) <- paste("day", now,
                     c("mutant", "control", "between"))

  bb <-  c(bb, aa)

  dist.list[[i]] <-  list(as.vector(between.dist),
                          as.vector(mu.dist),
                          as.vector(con.dist))
  names(  dist.list[[i]]) <- c("between", "mutant", "control")

  ##### ALPHA DIVERSITY

  alpha.t <- t.test(map.pick$invsimp ~ map.pick$genotype)
  asd.summary$alpha.t.pval[i] <- alpha.t$p.value
}

for(i in 1 : 5){
  # add bc to centroid to maps of each time point
  now <- time.line[i]

  otu.pick <- otu[map$time == now , ]
  map.pick <- map[map$time == now , ]

  assign(paste0("otu.", now), otu.pick)
  assign(paste0("map.", now), map.pick)
}

{

  par(mfrow=c(1, 1))
  par(mar=c(5.1, 5.1, 5.1, 5.1))

  p <- asd.summary$`PERMANOVA p val`

  plot(asd.summary$`PERMANOVA p val` ~ time.line, col = 1, lwd = 2,
       xlab = "Days post-weaning (DPW)",
       ylab = "", type = "l",   yaxt='n')

  points(p ~ time.line,
         pch = 21,
         bg = 1,
         cex = 2)

  abline (h = 0.05, lty = 2, col = 2)

  axis(2, at= p ,labels=round(p ,digits=2),
       col.axis=1, las=2, cex.axis=1, tck=-.01,
       col = 1)

  mtext("PERMANOVA p value", side=2, line=3, cex.lab=1,
        las=3, col=1)
  text(45, .055, "p value = 0.05", col = 2)

  par(new=TRUE)

  f <- asd.summary$`F stat`
  plot(asd.summary$`F stat` ~ time.line, type = "l",  xlab = "",
       ylab = "",  yaxt='n', col = 8, lwd = 2, lty = 2)

  points(asd.summary$`F stat` ~ time.line,
         pch = 21, bg = 8, cex = 2)

  axis(4, at=f,labels=round(f,digits=2),
       col.axis=1, las=2, cex.axis=1, tck=-.01, col = 1)
  mtext("PERMANOVA F statistic", side=4, line=3, cex.lab=1,
        las=3,  col = 1)

  legend("top",
         pch = 21,
         pt.cex = 2,
         pt.bg = c(1, 8),
         legend = c("p value", "F statistic"))
}
