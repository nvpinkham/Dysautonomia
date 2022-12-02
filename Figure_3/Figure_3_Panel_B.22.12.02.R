
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




asd.summary <- matrix(nrow = 5, ncol = 4)
colnames(asd.summary) <- c("PERMANOVA p val", "F stat", "DF", "R2")
rownames(asd.summary) <- sort(unique(map$time))
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
  asd.summary $DF[i] <- res$Df[1]
  asd.summary $R2[i] <- res$R2[1]


  asd.summary$n.F.wt[i] <- length(map.pick$sex[map.pick$sex == "f" & map.pick$genotype == "control"])
  asd.summary$n.M.wt[i] <- length(map.pick$sex[map.pick$sex == "m" & map.pick$genotype == "control"])

  asd.summary$n.F.mutants[i] <-length(map.pick$sex[map.pick$sex == "f" & map.pick$genotype == "mutant"])
  asd.summary$n.M.mutants[i] <- length(map.pick$sex[map.pick$sex == "m" & map.pick$genotype == "mutant"])
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


source_data_3b <- cbind(time.line,
                        asd.summary$`PERMANOVA p val`,
                        asd.summary$`F stat`)

colnames(source_data_3b) <- c("DPW", "PERMANOVA p val", "PERMANOVA F stat")

source_data_3b <- cbind("", "", source_data_3b)

source_data_3b[1,1] <- "Figure 3b"

write.csv(source_data_3b,"source_data/source_data_3b.csv", row.names = F)

write.table(asd.summary, "Statistical_summaries/Fig_3ab_tests.txt")


