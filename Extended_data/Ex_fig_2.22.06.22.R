

source("R/human_functions.22.06.18.R")
# run time : 10-20 minutes

map <- read.csv("data/Map_human_microbiome.csv")
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


f.res <- fisher.test(table(map$Disease.state, map$pam))
f.res

t.res <- t.test(map$invsimp ~ map$pam)
t.res

#Reason we did not use Chi squared: (25% of cells have < 5)
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5426219/#:~:text=The%20chi%2Dsquared%20test%20applies,especially%20for%20small%2Dsized%20samples.

map$dummy <- map$pam
map$dummy  <- map$dummy + sample(-100 : 100, nrow(map), replace = T)/2000

tab <- table(map$Disease.state, map$pam)


par(mfrow=c(2,1))

barplot(tab,
        col = c("#E69F00" ,"#56B4E9"),
        ylab = "sample count")

legend("topleft",
       fill = c( "#E69F00", "#56B4E9"),
       bty = "n",
       legend = c("Patient",
                  "Relative"))


title(paste0("Fisher's Exact Test\np-value =", round(f.res$p.value, 5),
             ", odds ratio = ", round(f.res$estimate, 5)))

vioplot(map$invsimp ~ map$pam, col = c(3, 4),
        drawRect = F,
        ylab = "Invsimp",
        xlab = "Cluster/Enterotype")


#points(map$invsimp ~ map$dummy, pch = 21, bg = map$col, cex = 1.5)
points(map$invsimp ~ map$dummy, pch = 21, bg = map$col)

legend("topleft",
       pch = 21,
       pt.bg = c( "#E69F00", "#56B4E9"),
       bty = "n",
       legend = c("Patient",
                  "Relative"))

title(paste0("T-test\np-value =", round(t.res$p.value, 10),
             ", effect size= ", round(diff(t.res$estimate), 5)))

