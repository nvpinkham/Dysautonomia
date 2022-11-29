getwd()
source("R/human_functions.22.11.27.R")

map <- read.csv("data/Map_human_microbiome.22.11.21.csv")
row.names(map) <- map$sample.id
otu <- read.csv("data/FD_OTU_human21.csv", row.names = 1)
tax <- read.csv("data/FD_OTU_human21_taxonomy.csv")

map$enterotype[map$pam == 1] <- "Cluster1"# more FD
map$enterotype[map$pam != 1] <- "Cluster2"# less FD

map$discription <- map$enterotype

map$Family.ID[map$pam == 1 & map$Disease.state == "Relative"]
map.1.r <- map[map$pam == 1 & map$Disease.state == "Relative", ]

boxplot(map$invsimp ~ map$pam)
res <- t.test(map$invsimp ~ map$pam)
diff(res$estimate)

par(mar=c(4.1, 4.1, 5.1, 1))

rf.res <- pairwiseRF(otu, map, dis_1 = "Cluster1", dis_2 = "Cluster2", Factor = 2.5, plot = T)

comps <- c("Cluster1", "Cluster2")

aa <- match(names(rf.res[[1]]), colnames(otu))
bb <- match(names(rf.res[[2]]), colnames(otu))

meta.dif <- c(aa, bb)

colnames(otu) <- paste(colnames(otu), "\n", tax$genus)

otu.counts <- otu

otu <- log(otu)
otu [otu == -Inf] <- 0

m1 <- otu[map$discription == comps[1],
          aa]

m1 <- cbind(m1, otu[map$discription == comps[1],
                    bb])


m2 <- otu[map$discription == comps[2],
          aa]

m2 <- cbind(m2, otu[map$discription ==  comps[2],
                    bb])

m1 <- m1[,ncol(m1):1]
m2 <- m2[,ncol(m2):1]

a <- 0 : (length(m1) -1) * 4

par(mar=c(5.1, 13.1, 5.1, 1))

#setEPS()
#postscript(paste0("Figure2_PanelA",  Sys.Date(),".eps"), width = 6, height = 6)

par(mar=c(4.1, 14.1, 5.1, 8.1))

boxplot(m1, las = 2, horizontal = TRUE,
        at = a,
        xlim = c(-2, max(a) + 6),
        pch = 21,
        cex = .5,
        ylim = range(c(m1, m2)),
        xlab = "\n natural log OTU count",
        bg = "gray40",
        col = "gray40")

boxplot(m2 , las = 2, horizontal = TRUE, yaxt="none",
        at = a  - 1, add= T, names = rep("", length(a)),
        pch = 21,
        cex = .5,
        bg = "gray80",
        col = "gray80")

abline(h = length(bb) * 4 - 2, col = "red", lty = 2, lwd = 2)
#abline(h = length(bb) * 4 - 2, col = "gray30", lty = 2, lwd = 2)

legend("top", fill= c("gray40", "gray80"),
       legend = c("Cluster/Enterotype 1", "Cluster/Enterotype 2"))

#dev.off()

otu.res <- as.data.frame(matrix(nrow = length(m1), ncol = 7))
colnames(otu.res) <- c("P_val", "df", "effect_size", "t.stat",
                       "conf_int_lo",  "conf_int_hi", "wil_p")

for(i in 1 : length(m1)){

  res <- t.test(m1[,i], m2[,i])

  otu.res$P_val[i] <- res$p.value
  otu.res$df[i] <- res$parameter
  otu.res$t.stat[i] <- res$statistic
  otu.res$effect_size[i] <- diff(res$estimate)
  otu.res$conf_int_lo[i] <- res$conf.int[1]
  otu.res$conf_int_hi[i] <- res$conf.int[2]
}


table(map$enterotype, map$Disease.state)

##################### compile source data file ###########################
source_data_2a <- cbind(map$sample.id, map$pam)
colnames(source_data_2a) <- c("sample.id", "enterotype/cluster")
source_data_2a <- cbind("", source_data_2a, "","", as.matrix(otu.counts))
source_data_2a[1,1] <- "Figure 2a"
source_data_2a[1,5] <- "OTU counts"

write.csv(source_data_2a,"source_data/source_data_2a.csv", row.names = F)












