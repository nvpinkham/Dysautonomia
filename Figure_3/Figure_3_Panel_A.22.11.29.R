
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

######################### Expand MAP file ############################
# mice were separated by sex and genotype
map$cage <- paste0(map$genotype, ".", map$sex)

map$invsimp <- diversity(otu, "invsimp")
m.dist <- vegdist(otu)
map$discription <- paste("day", map$time, map$genotype)

#map$color[map$genotype == "mutant"] <- "#E69F00" # colors to match human data
#map$color[map$genotype != "mutant"] <- "#56B4E9"

map$color[map$genotype == "mutant"] <-  "#D55E00"
map$color[map$genotype != "mutant"] <-  "#F0E442"

par(mfrow=c(1, 3))
par(mar=c(4.1, 4.1, 4.1, 1))

time.pick <- c(0, 25, 54)

for(i in 1 : length(time.pick)){

  now <- time.pick[i]

  otu.pick <- otu[map$time == now , ]
  map.pick <- map[map$time == now , ]

  set.seed(42)
  pick.bc <- vegdist(otu.pick)
  otu.nmds <- metaMDS(pick.bc, plot = F,  k = 2, try = 20)

  {plot(otu.nmds$points,
        main = paste("Day", time.pick[i],
                     "\nn =", nrow(otu.pick)))

    ordispider(otu.nmds$points, groups = map.pick$discription,
               col = 8, lty = 1)
    points(otu.nmds$points, pch = 21, bg = map.pick$color, cex = 2)

    #### compile sorce data
    source_data_3a <- cbind(map.pick$sample.id,
                            map.pick$genotype)

    colnames(source_data_3a) <- c("sample.id", "Genotype")

    source_data_3a <- cbind("", "", source_data_3a, "", "", otu.nmds$points)
    source_data_3a[1,2] <- paste("Day", now)
    source_data_3a[1,6] <- "NMDS points"

    source_data_3a[1,8] <- "Bray-Curtis Dissimilarity between samples"

    pick.bc <- as.matrix(pick.bc)
    pick.bc[upper.tri(pick.bc, diag = T)] <- ""
    source_data_3a <- cbind(source_data_3a, map.pick$sample.id, pick.bc)

    source_data_3a <- rbind(source_data_3a, "", "")

    ###################

    if(now == 0){

      legend("bottomleft",
             pch = 21,
             cex = 1.5,
             pt.bg = c(  map.pick[1,]$color,
                         map.pick[2,]$color),
             legend = c(  map.pick[1,]$genotype,
                          map.pick[2,]$genotype))
      source_data_3a[1,1] <- "Figure 3a"

    }

    write.csv(source_data_3a,
              paste0("source_data_3_dummy.", now,".csv"),
              row.names = F)

  }
}


################################################################################

system("cat source_data_3_dummy* > Statistical_summaries/source_data_3a.csv")
system("rm source_data_3_dummy*")
