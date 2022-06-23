getwd()
source("R/human_functions.22.06.18.R")

map <- read.csv("data/Map_human_microbiome.csv")
row.names(map) <- map$sample.id
otu <- read.csv("data/FD_OTU_human21.csv", row.names = 1)
tax <- read.csv("data/FD_OTU_human21_taxonomy.csv")

map$enterotype[map$pam == 1] <- "FD"
map$enterotype[map$pam != 1] <- "Normal"

map$discription <- map$enterotype

map$Family.ID[map$pam == 1 & map$Disease.state == "Relative"]
map.1.r <- map[map$pam == 1 & map$Disease.state == "Relative", ]

t.test(map$invsimp ~ map$pam)

par(mar=c(4.1, 4.1, 5.1, 1))

rf.res <- pairwiseRF(otu, map, dis_1 = "FD", dis_2 = "Normal", Factor = 2.5, plot = T)


comps <- c("FD", "Normal")

aa <- match(names(rf.res[[1]]), colnames(otu))
bb <- match(names(rf.res[[2]]), colnames(otu))

meta.dif <- c(aa, bb)

colnames(otu) <- paste(colnames(otu), "\n", tax$genus)

otu <- log(otu)
otu [otu == -Inf] <- 0

m1 <- otu[map$discription == comps[1],
          aa]

m1 <- cbind(m1, otu[map$discription == comps[1],
                    bb])


m2 <- otu[map$discription ==  comps[2],
          aa]

m2 <- cbind(m2, otu[map$discription ==  comps[2],
                    bb])

m1 <- m1[,ncol(m1):1]
m2 <- m2[,ncol(m2):1]

a <- 0 : (length(m1) -1) * 4

par(mar=c(5.1, 13.1, 5.1, 1))

setEPS()       # Set postscript arguments
postscript(paste0("Figure2_PanelA",
                  #  Sys.Date(),
                  ".eps"),
           width = 6, height = 6)

par(mar=c(4.1, 14.1, 5.1, 1))

boxplot(m1, las = 2, horizontal = TRUE,
        at = a,
        xlim = c(-2, max(a) + 6),
        pch = 21,
        cex = .5,
        bg = 3,
        ylim = range(c(m1, m2)),
        xlab = "\n natural log OTU count",
        col = 3)

boxplot(m2 , las = 2, horizontal = TRUE, yaxt="none",
        at = a  - 1, add= T, names = rep("", length(a)),
        pch = 21,
        cex = .5,
        bg = "dodgerblue",
        col = "dodgerblue")

abline(h = length(bb) * 4 - 2, col = 2)


legend("top", fill= c(3, "dodgerblue"),
       legend = c("Enterotype 1", "Enterotpe 2"))

dev.off()






