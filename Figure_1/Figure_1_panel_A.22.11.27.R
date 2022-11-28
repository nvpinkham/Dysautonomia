
source("R/human_functions.22.11.27.R")

# make Ordination of beta diversity in human participants
# should take under a minute

# import data
map <- read.csv("data/Map_human_microbiome.22.11.21.csv")
row.names(map) <- map$sample.id
otu <- read.csv("data/FD_OTU_human21.csv", row.names = 1)

table(map$Disease.state)
length(unique(map$Family.ID))

reconcile(map, otu) # make sure samples match

set.seed(42)
otu.nmds <- metaMDS(otu, distance = "bray") # make ordination

otu.bc <- vegdist(otu) # calc dissimilarity

set.seed(42)
res <- adonis2(otu.bc ~ map$Disease.state)# Run permanova compairing patients to relatives
# pairing not incorperated, some subjects resampled

res1 <- paste("BC dist ~ disease state \n f stat=",
              round(res$F[1], 3),
              "p val=", res$`Pr(>F)`[1],
              "\ndf =", res$Df[1],
              "R2 = ",  round(res$R2[1], 5))

plot(otu.nmds$points, cex = .2, main = res1)

num.pairs <- NULL

for(i in unique(map$Family.ID)){
  # draw lines connecting all samples to their possible paired samples (samples from opposite category in the same family)

  nmds.pick <- otu.nmds$points[map$Family.ID == i,]
  map.pick <- map[map$Family.ID == i,]

  if(length(unique(map.pick$Disease.state)) > 1){

    nmds.pick.sick <- matrix(nmds.pick[map.pick$Disease.state == "Patient",],
                                ncol = 2, byrow = F)
    nmds.pick.control <- matrix(nmds.pick[map.pick$Disease.state == "Relative",],
                                ncol = 2, byrow = F)

    for(j in 1 : nrow(nmds.pick.sick)){

      x0 <- nmds.pick.sick[j, 1]
      y0 <- nmds.pick.sick[j, 2]

      for(k in 1: nrow(nmds.pick.control)){

        arrows(x0,
               y0,
               nmds.pick.control[k,1],
               nmds.pick.control[k,2],
               col = 8, angle = 0)

        num.pairs <- c(num.pairs, 1)
      }
    }
  }
}

length(num.pairs)

points(otu.nmds$points,
       cex =  map$invsimp / 10 + .2,
       bg = map$col,
       pch = 21)

legend("bottomright",
       pt.bg = map$col,
       pch = 21, legend =
         c("patient", "relative"))

legend(pch = 21, pt.bg  = 8,
       pt.cex = c(1, 3),
       legend = c("low alpha", "high alpha"), "topright")

a <- lme4::lmer(map$invsimp ~ map$Disease.state + map$Family.ID +  (1|map$DonorID))
b <- lmer(map$invsimp ~ map$Disease.state + (1|map$DonorID))
#b <- lm(map$invsimp ~ )
anova(a, b)

summary(a)
confint(a, level = 0.95)

b <- lmer(map$invsimp ~  map$Family.ID + (1|map$DonorID))
anova(a, b)

################################################################################

source_data.1a <- cbind(map$sample.id,  map$DonorID, map$Family.ID, map$Disease.state, map$invsimp)

colnames(source_data.1a) <- c("sample.id",  "DonorID", "Family.ID", "Disease.state", "Invsimp")

source_data.1a <- cbind("", source_data.1a, "", "")
source_data.1a[1,1] <- "Figure 1a"
source_data.1a[1,8] <- "Bray-Curtis Dissimilarity between samples"

otu.bc <- as.matrix(otu.bc)
otu.bc[upper.tri(otu.bc, diag = T)] <- ""
source_data.1a <- cbind(source_data.1a, map$sample.id, otu.bc)

write.csv(source_data.1a,"source_data_1a.csv", row.names = F)

res

