
source("R/human_functions.22.11.27.R")
# run time : 10-20 minutes

map <- read.csv("data/Map_human_microbiome.22.11.21.csv")
row.names(map) <- map$sample.id
otu <- read.csv("data/FD_OTU_human21.csv", row.names = 1)

table(map$Disease.state)
length(unique(map$Family.ID))

reconcile(map, otu)

sils <- NULL
sils.individual <- matrix(ncol = 15, nrow = nrow(otu))
sils.group <- matrix(ncol = 15, nrow = 2)

cw <- list()

otu.bc <- vegdist(otu)

for(i in 2:15){
  #calc sil width for whole dataset and by group

  pam.i <- pam(otu.bc, k = i)
  sils[i] <- pam.i$silinfo$avg.width

  a <- match( map$sample.id, row.names( pam.i$silinfo$widths))
  b <- pam.i$silinfo$widths[a,]

  p <- mean(b[map$Disease.state == "Patient" , 3])
  r <- mean(b[map$Disease.state != "Patient" , 3])

  sils.individual[,i] <- b[,3]

  sils.group[1,i] <- p
  sils.group[2,i] <- r

  cw[[i]] <- pam.i$silinfo$clus.avg.widths
}


plot(sils, type = "l", lty = 1, col = 8, lwd = 2,
     xlim = c(2,15),
     # ylim = range(sils.individual, na.rm = T),
     ylim = range(sils.group, na.rm = T),
     #ylim = range(cw, na.rm = T),

     xlab = "Number of Clusters (k)",
     ylab = "Average silhouette width")
points(sils, pch =19, col = 8, cex = 2)

points(sils.group[1,], type = "l", lwd = 2, col = "#E69F00")
# this is not entrotype,
# this is the avege of all patient's sil widths
#  -REGARDLESS of enterotype/Pam cluster
points(sils.group[1,], pch = 19, col ="#E69F00" , cex = 2)

points(sils.group[2,], type = "l",lwd = 2, col = "#56B4E9")
points(sils.group[2,], pch = 19, col = "#56B4E9", cex = 2)

#for(i in 2 : length(cw)){
#  points(x= rep(i, i), y = cw[[i]], pch = "+" )
#}

legend("topright",
       col = c( "#E69F00", "#56B4E9", 8),
       pch = 19,
       bty = "n",
       legend = c("Average Patient",
                  "Average Relative",
                  "Average All"))

pam.2 <- pam(otu.bc, k =2)
sil.2 <- silhouette(pam.2)
a <- match(row.names(pam.2$silinfo$widths), map$sample.id)
map$sample.id[a]

plot(sil.2, col = map$col[a])
legend("bottom",
       fill = c( "#E69F00", "#56B4E9"),
       bty = "n",
       legend = c("Patient",
                  "Relative"))

res <- fisher.test(table(map$pam, map$Disease.state), alternative = "two.sided")

####################### EXPORT SORCE DATA FILE #################################

sup_data_2a <- cbind(1:15, sils, sils.group[1,], sils.group[2,])[-1,]

colnames(sup_data_2a) <- c("K", "Average sil width",
                           "Average sil width only patients",
                           "Average sil width only relatives")

sup_data_2a <- cbind("", "", sup_data_2a)
sup_data_2a[1,1] <- "Supplemental figure 2a"

write.csv(sup_data_2a,"sup_data_2a.csv", row.names = F)

sup_data_2b <- cbind(map$sample.id[a], map$Disease.state[a], sil.2)

colnames(sup_data_2b)[2] <- c("Disease State")

sup_data_2b <- cbind("", "", sup_data_2b)
sup_data_2b[1,1] <- "Supplemental figure 2b"

write.csv(sup_data_2b,"source_data/sup_data_2b.csv", row.names = F)

res <- capture.output(res)
writeLines(res, "Statistical_summaries/Sup_fig_2.txt")



