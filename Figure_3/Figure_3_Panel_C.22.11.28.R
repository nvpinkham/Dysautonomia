
library(vegan)
library(vioplot)

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
otu.dist <- vegdist(otu)
b.res <- betadisper( otu.dist, paste(map$time, map$genotype))

map$discription
map$discription <- paste("  Day", map$time, "\n", map$genotype)

# setEPS()       # Set postscript arguments
# postscript(paste0("Figure3_panelC_", Sys.Date(), ".eps"),
#            width = 10, height = 6)

vioplot(b.res$distances ~  map$discription,
        drawRect = F,
        font = 2,
        col= c("#F0E442", "#D55E00"),
        xlab = "Day Post Weaning (DPW)",
        ylab = "Bray-Curtis Distance to Group Centroid")

abline(v = c(2.5, 4.5, 6.5, 8.5), lty = 2)


tps <- sort(unique(map$time))
tps

t.res <- list()

for(i in 1:5){

  t.res[[i]] <- t.test(b.res$distances[map$time == tps[i]] ~  map$discription[map$time == tps[i]])

}

segments(9, .36, 10, .36)
text(9.5, .37, "p < 0.05", font = 4)

title.t <- paste("Day 54 t-test\neffect size =",
                 round(diff(t.res[[5]]$estimate), 5),
                 "p val = ", round(t.res[[5]]$p.value, 5),
                 "\n95% CI =", round(t.res[[5]]$conf.int[1], 5),
                 ",", round(t.res[[5]]$conf.int[2], 5),
                 "\ndf = ", round(t.res[[5]]$parameter, 5))


title(title.t)

agg <- aggregate(b.res$distances, list(map$discription), mean)
agg

for(i in 1 : 10){
  segments(i-.3, agg$x[i], i + .3, agg$x[i], lwd = 2)
}

mice <- unique(map$mouse)

days <- c(0, 15, 25, 40, 54)
con.x <- seq(1, 9, 2)
mu.x <- seq(2, 10, 2)

for(i in 1 : length(mice)){

  y <- b.res$distances[map$mouse == mice[i]]

  x <- map$time[map$mouse == mice[i]]

  if(map$genotype[map$mouse == mice[i]][1] == "control"){

    x <- con.x[which(days %in% x)]
    lines(x, y, col = 8 , lty = 2)

  }else{

    x <- mu.x[which(days %in% x)]
    lines(x, y, col = 8, lty = 4)

  }
}


groups <- sort(unique(map$discription))

for(i in 1 : length(groups)){

  y <- b.res$distances[map$discription == groups[i]]

  points(rep(i, length(y)) , y, pch = 21, bg = 8)
}

# dev.off()

all(names(b.res$distances) == map$sample.id)
source_data_3c <- cbind(map$sample.id, map$mouse, map$genotype, map$time, b.res$distances)

colnames(source_data_3c) <- c("sample.id",  "Mouse", "genotype", "time point", "Bray-Curtis Distance to Group Centroid")

source_data_3c <- cbind("", "", source_data_3c)
source_data_3c[1,1] <- "Figure 3c"

write.csv(source_data_3c,"source_data_3c.csv", row.names = F)
t.res[[5]]

