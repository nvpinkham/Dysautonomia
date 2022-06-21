
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

######################### FINISH MAP ############################
# mice were separated by sex and genotype
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





par(mfrow=c(1, 3))
par(mar=c(4.1, 4.1, 4.1, 1))

time.pick <- c(0, 25, 54)

for(i in 1 : length(time.pick)){

  now <- time.pick[i]

  otu.pick <- otu[map$time == now , ]
  map.pick <- map[map$time == now , ]

  set.seed(42)
  otu.nmds <- metaMDS(otu.pick, plot = F,  k = 2, try = 20)

  {plot(otu.nmds$points,
        main = paste("Day", time.pick[i],
                     "\nn =", nrow(otu.pick)))

    ordispider(otu.nmds$points, groups = map.pick$discription,
               col = 8, lty = 1)
    points(otu.nmds$points, pch = 21, bg = map.pick$color, cex = 2)

    if(now == 0){
      legend("bottomleft",
             pch = 21,
             cex = 1.5,
             pt.bg = c(  map.pick[1,]$color,
                         map.pick[2,]$color),
             legend = c(  map.pick[1,]$genotype,
                          map.pick[2,]$genotype))
    }
  }
}

