
source("R/human_functions.22.11.21.R")
setwd("/Users/nickpinkham/Desktop/Dysautonomia")

meta.serum <- read.csv("data/HumanSerumUnpaired_biomassnorm_metalabels.csv", row.names = 1)
map.serum <- read.csv("data/Map_human_metabolome_serum.csv", row.names = 1)

serum.norm <- norm.human(metabolites = meta.serum, samples.pick = row.names(meta.serum))

serum.norm$Label <- NULL
s.d <- dist(serum.norm)

serum.norm$X4.Pyridoxate <- NULL # this was not a real finding

serum.plsda <- mixOmics::splsda(serum.norm, map.serum$Disease.state, ncomp = 3)

sigs <- round(serum.plsda$prop_expl_var$X * 100)
sigs <- paste0("component ", 1:3, ": ", sigs, "%")

map.serum.plsda <- with(map.serum, vegan3d::ordiplot3d(serum.plsda$variates$X,
                                                       pch = 21,
                                                       ax.col = NULL,
                                                       bg = map.serum$col,
                                                       xlab = sigs[1],
                                                       ylab = sigs[2],
                                                       zlab = sigs[3]))



# setEPS()       # Set postscript arguments
# postscript(paste0("Serum_metabolite_plsda_.",Sys.Date(), ".eps"),width = 5, height = 5)
scatterplot3d(serum.plsda$variates$X,
              pch = 21,
              ax.col = NULL,
              bg = map.serum$col,
              xlab = sigs[1],
              ylab = sigs[2],
              zlab = sigs[3])

loadings <- map.serum.plsda$points
#loadings <- serum.plsda$variates$X

for(i in map.serum$Family.ID){


  loadings.pick <- loadings[map.serum$Family.ID == i,]
  map.pick <- map.serum[map.serum$Family.ID == i,]

  if(length(unique(map.pick$Disease.state)) > 1){


    loadings.pick.sick <- matrix(loadings.pick[map.pick$Disease.state == "Patient",], ncol = 2, byrow = F)
    loadings.pick.control <- matrix(loadings.pick[map.pick$Disease.state == "Relative",],ncol = 2, byrow = F)

    for(j in 1 : nrow(loadings.pick.sick)){

      x0 <- loadings.pick.sick[j, 1]
      y0 <- loadings.pick.sick[j, 2]

      for(k in 1: nrow(loadings.pick.control)){
        arrows(x0,
               y0,
               loadings.pick.control[k,1],
               loadings.pick.control[k,2],
               col = 8, angle = 0)


      }
    }
  }
}

points(loadings,
       pch = 21,
       bg = map.serum$col)

set.seed(42)
res <- adonis2(dist(serum.norm) ~ map.serum$Disease.state)
res <- paste("distance ~ disease state \n f stat",
             round(res$F[1], 3), "p val", res$`Pr(>F)`[1])
title(res)

#dev.off()

source_data.1f <- cbind(map.serum$sample.id,  map.serum$DonorID, map.serum$Family.ID, map.serum$Disease.state)

colnames(source_data.1f) <- c("sample.id",  "DonorID", "Family.ID", "Disease.state")

source_data.1f <- cbind("", source_data.1f, "", "", serum.plsda$variates$X)
source_data.1f[1,1] <- "Figure 1e"
source_data.1f[1,7] <- "PLSDA coordinates"

write.csv(source_data.1f,"source_data_1f.csv", row.names = F)
