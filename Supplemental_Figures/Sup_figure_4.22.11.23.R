source("R/human_functions.22.11.21.R")
getwd()

map <- read.csv("data/FDMS_map_human22.csv", row.names = 1)
otu <- read.csv("data/FDMS_OTU_human22.csv", row.names = 1)
#tax <- read.csv("data/FDMS_OTU_human22_taxonomy.csv", row.names = 1)

map.isme <- read.csv("data/isme19.csv")
map.fd <- read.csv("data/Map_human_microbiome.22.11.21.csv")
map.fd <- map.fd[map.fd$sample.id != "FD006" , ]

map$donor <- map$sample.id
map$donor[match(map.isme$Experiment.Title, map$sample.id)] <- map.isme$donor

map$donor[match(map.fd$sample.id, map$sample.id)] <- map.fd$DonorID

map$discription <- map$Disease.state
map$id <- map$Experiment.Title
all(map$id == row.names(otu))

map$discription[map$Study.Title == "Temporal dynamics of Enterobacteriaceaein healthy human adults"] <- "Martinson"
map$discription[map$Study.Title == "Microbiome analysis complements fecal occult blood test for improved detection"] <- "Baxter"

map$col[map$discription == "Martinson"] <- "tomato"
map$col[map$discription == "Baxter"] <- "green"
map$col[map$Disease.state == "Relative"] <- "#56B4E9"
map$col[map$Disease.state == "Patient"] <- "#E69F00"

map$invsimp <- diversity(otu, "invsimp")

res2title <- function(res){
  res <- paste("BC dist ~ group p val", res$`Pr(>F)`[1])
  res
}


groups <- c("Patient", "Martinson")

make_nmds_FDMS <- function(groups){

  pick <- map$discription %in% groups

  res <- adonis2(otu[pick,] ~ map$discription[pick])# Run

  set.seed(42)
  otu.nmds <- metaMDS(otu, distance = "bray") # make ordination

  ord.dist <- dist(otu.nmds$points)
  bc.dist <- vegdist(otu)


  plot(otu.nmds$points,
       xlim = c(min(otu.nmds$points[,1]),
                max(otu.nmds$points[,1])) * 1.5, cex = map$invsimp / 10 + .2, col = 8)
  # main = res)

  title(paste(toString(groups), res2title(res)))

  map$invsimp <- diversity(otu, "invsimp")

  points(otu.nmds$points[pick,],
         cex =  map$invsimp[pick] / 10 + .2,
         bg = map$col[pick],
         pch = 21)

  leg1 <- aggregate(map$col, list(map$discription), unique)

  legend("bottomleft",
         pt.bg = leg1$x,
         pch = 21, legend = leg1$Group.1)

  legend(pch = 21, pt.bg  = 8,
         pt.cex = c(1, 3),
         legend = c("low alpha", "high alpha"), "topleft")

  return(otu.nmds$points)
}



map$discription2 <- map$discription
map$col2 <- map$col

map$discription <- map$Disease.state
map$col[map$Disease.state == "Healthy"] <- "purple"


otu.nmds <- make_nmds_FDMS(groups = c("Healthy", "Relative"))

sup_data_4 <- cbind(map$sample.id, map$Sample.Accession, map$discription, map$discription2)
sup_data_4[!sup_data_4[,4] %in% c("Martinson", "Baxter") , 4] <- "Pinkham"
colnames(sup_data_4) <- c("sample.id", "sample.accession", "group", "study")

sup_data_4 <- cbind("", "", sup_data_4,"", "", otu.nmds, "", "", as.matrix(otu))
sup_data_4[1,1] <- "Supplemental figure 4a & 4b"
sup_data_4[1,2] <- "sample information"
sup_data_4[1,8] <- "NMDS points"
sup_data_4[1,12] <- "OTU counts"

write.csv(sup_data_4,"sup_data_4ab.csv", row.names = F)
################################################################################

make_nmds_FDMS(groups = c("Healthy", "Patient"))

################################################################################

map$discription <- map$discription2
map$col <- map$col2

groups <- unique(map$discription)

for(i in 1 : 4){
  otu.pick <- otu[map$discription == groups[i],]
  assign(paste0("centroid.", groups[i]), colMeans(otu.pick))
}

leg1 <- aggregate(map$col, list(map$discription), unique)
leg1$Group.1

cents <- rbind(centroid.Baxter,
               centroid.Martinson,
               centroid.Patient,
               centroid.Relative)

row.names(cents) <- gsub("centroid.", "", row.names(cents))

for(i in 1 : 4){
  otu.pick <- otu[map$discription == groups[i],]
  assign(paste0("centroid.", groups[i]), colMeans(otu.pick))
}

otu.bc <- vegdist(otu)

beta <- betadisper(otu.bc, map$discription)
beta$centroids

dists <- vegdist(cents)

par(mfrow=c(2,1))

par(xpd = TRUE) #Draw outside plot area
library(labdsv)
pco <- pco(dists)

plot(pco, xlim = c(-.5, .3), pch = 21, bg =  leg1$x, cex = 2)
text(pco$points, row.names(pco$points), pos = 3)

dists2 <- spaa::dist2list(dists)

dists2 <- dists2[!duplicated(dists2$value),]
dists2 <- dists2[dists2$value != 0,]

dists2$p1.y <- dists2$p1.x <- NA
dists2$p2.y <- dists2$p2.x <- NA

for(i in 1:4){
  dists2$p1.x[dists2$col == row.names(pco$points)[i]] <- pco$points[i,1]
  dists2$p1.y[dists2$col == row.names(pco$points)[i]] <- pco$points[i,2]

  dists2$p2.x[dists2$row == row.names(pco$points)[i]] <- pco$points[i,1]
  dists2$p2.y[dists2$row == row.names(pco$points)[i]] <- pco$points[i,2]

}

segments(x0 = dists2$p1.x, y0 = dists2$p1.y, x1 = dists2$p2.x, y1 = dists2$p2.y, col = 8)


for(i in 1 : 6){
  #segments(x0 = dists2$p1.x[i], y0 = dists2$p1.y[i], x1 = dists2$p2.x[i], y1 = dists2$p2.y[i], col = 2)

  dists2[i,]
  points(mean(c(dists2$p1.x[i], dists2$p2.x[i])),
         mean(c(dists2$p1.y[i], dists2$p2.y[i])), col = 8, pch = 19)

  text(mean(c(dists2$p1.x[i], dists2$p2.x[i])),
       mean(c(dists2$p1.y[i], dists2$p2.y[i])),
       round(dists2$value[i], 4), pos = 4, col = 8)
}


################################################################################
dists <- as.matrix(dists)
dists[upper.tri(dists, diag = T)] <- ""

group <- row.names(dists)
sup_data_4d <- cbind("", "", group, dists, "", "", group, as.matrix(cents))
sup_data_4d[1,1] <- "Supplemental figure 4d"
sup_data_4d[1,2] <- "Bray-Curtis distance between group centriods"
sup_data_4d[1,9] <- "Mean OTU counts by group"

write.csv(sup_data_4d,"sup_data_4d.csv", row.names = F)

################################################################################

groups <- unique(map$discription)

for(i in 1 : 4){
  otu.pick <- otu[map$discription == groups[i],]

  Baxter <- as.matrix(vegdist(rbind(centroid.Baxter, otu.pick)))[-1, 1]
  Martinson <- as.matrix(vegdist(rbind(centroid.Martinson, otu.pick)))[-1, 1]
  Relatives <- as.matrix(vegdist(rbind(centroid.Relative, otu.pick)))[-1, 1]
  Patients <- as.matrix(vegdist(rbind(centroid.Patient, otu.pick)))[-1, 1]

  c(Baxter, Martinson, Relatives, Patients)

  # boxplot(Baxter, Martinson, Relatives, Patients, main = groups[i],         names = c("Baxter", "Martinson", "Relatives", "Patients"), ylim = c(0, 1))



  t.test(Baxter, Martinson, paired = T)
  t.test(Baxter, Relatives, paired = T)
  t.test(Baxter, Patients, paired = T)

  t.test(Martinson, Relatives, paired = T)
  t.test(Martinson, Patients, paired = T)

  t.test(Martinson, Patients, paired = T)

}

leg1$Group.1

cents <- rbind(centroid.Baxter,
               centroid.Martinson,
               centroid.Patient,
               centroid.Relative)

row.names(cents) <- gsub("centroid.", "", row.names(cents))

#####################

groups <- unique(map$Disease.state)

for(i in 1 : 3){
  otu.pick <- otu[map$Disease.state == groups[i],]
  assign(paste0("centroid.", groups[i]), colMeans(otu.pick))
}


cents <- rbind(centroid.Healthy,
               centroid.Patient,
               centroid.Relative)

row.names(cents) <- gsub("centroid.", "", row.names(cents))

for(i in 1 : 3){
  otu.pick <- otu[map$Disease.state == groups[i],]
  assign(paste0("centroid.", groups[i]), colMeans(otu.pick))
}

dists <- vegdist(cents)

library(labdsv)
pco <- pco(dists)

#pco$points[,2] <- pco$points[,2] * -1
#pco$points[,1] <- pco$points[,1] * -1

plot(pco, xlim = c(-.5, .3), pch = 21, bg =  c("purple", "#E69F00", "#56B4E9"), cex = 2)
text(pco$points, row.names(pco$points), pos = 3)


dists2 <- spaa::dist2list(dists)

dists2 <- dists2[!duplicated(dists2$value),]
dists2 <- dists2[dists2$value != 0,]

dists2$p1.y <- dists2$p1.x <- NA
dists2$p2.y <- dists2$p2.x <- NA

for(i in 1:3){
  dists2$p1.x[dists2$col == row.names(pco$points)[i]] <- pco$points[i,1]
  dists2$p1.y[dists2$col == row.names(pco$points)[i]] <- pco$points[i,2]

  dists2$p2.x[dists2$row == row.names(pco$points)[i]] <- pco$points[i,1]
  dists2$p2.y[dists2$row == row.names(pco$points)[i]] <- pco$points[i,2]

}


segments(x0 = dists2$p1.x, y0 = dists2$p1.y, x1 = dists2$p2.x, y1 = dists2$p2.y, col = 8)


for(i in 1 : 3){
  #segments(x0 = dists2$p1.x[i], y0 = dists2$p1.y[i], x1 = dists2$p2.x[i], y1 = dists2$p2.y[i], col = 2)

  dists2[i,]
  points(mean(c(dists2$p1.x[i], dists2$p2.x[i])),
         mean(c(dists2$p1.y[i], dists2$p2.y[i])), col = 8, pch = 19)

  text(mean(c(dists2$p1.x[i], dists2$p2.x[i])),
       mean(c(dists2$p1.y[i], dists2$p2.y[i])),
       round(dists2$value[i], 4), pos = 4, col = 8)
}

################################################################################
dists <- as.matrix(dists)
dists[upper.tri(dists, diag = T)] <- ""

group <- row.names(dists)
sup_data_4c <- cbind("", "", group, dists, "", "", group, as.matrix(cents))
sup_data_4c[1,1] <- "Supplemental figure 4d"
sup_data_4c[1,2] <- "Bray-Curtis distance between group centriods"
sup_data_4c[1,8] <- "Mean OTU counts by group"

write.csv(sup_data_4c,"sup_data_4c.csv", row.names = F)
################################################################################
