
source("R/human_functions.22.11.27.R")

meta.stool <- read.csv("data/HumanStoolUnpaired_biomassnorm_metalabels.csv", row.names = 1)
map.stool <- read.csv("data/Map_human_metabolome_stool.22.11.21.csv", row.names = 1)

stool.norm <- norm.human(metabolites = meta.stool, samples.pick = row.names(meta.stool))

stool.norm$Label <- NULL

stool.plsda <- mixOmics::splsda(stool.norm, map.stool$Disease.state, ncomp = 3)

sigs <- round(stool.plsda$prop_expl_var$X * 100)
sigs <- paste0("component ", 1:3, ": ", sigs, "%")

map.stool.plsda <- with(map.stool, vegan3d::ordiplot3d(stool.plsda$variates$X,
                                                       pch = 21,
                                                       ax.col = NULL,
                                                       bg = map.stool$col,
                                                       xlab = sigs[1],
                                                       ylab = sigs[2],
                                                       zlab = sigs[3]))


#setEPS()
#postscript(paste0("Stool_metabolite_plsda_.", Sys.Date(),".eps"),width = 5, height = 5)

scatterplot3d(stool.plsda$variates$X,
              pch = 21,
              ax.col = NULL,
              bg = map.stool$col,
              xlab = sigs[1],
              ylab = sigs[2],
              zlab = sigs[3])

loadings <- map.stool.plsda$points

num.pair <- NULL

for(i in map.stool$Family.ID){


  loadings.pick <- loadings[map.stool$Family.ID == i,]
  map.pick <- map.stool[map.stool$Family.ID == i,]

  if(length(unique(map.pick$Disease.state)) > 1){


    loadings.pick.sick <- matrix(loadings.pick[map.pick$Disease.state == "Patient",], ncol = 2, byrow = F)
    loadings.pick.control <- matrix(loadings.pick[map.pick$Disease.state == "Relative",],ncol = 2, byrow = F)

    for(j in 1 : nrow(loadings.pick.sick)){

      x0 <- loadings.pick.sick[j, 1]
      y0 <- loadings.pick.sick[j, 2]

      for(k in 1: nrow(loadings.pick.control)){

        num.pair <- c(num.pair, 1)
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
       bg = map.stool$col)

set.seed(42)
res <- adonis2(dist(stool.norm) ~ map.stool$Disease.state)
res1 <- paste("distance ~ disease state \n f stat=",
              round(res$F[1], 3),
              "p val=", res$`Pr(>F)`[1],
              "\ndf =", res$Df[1], "R2 = ",  round(res$R2[1], 5))
title(res1)

map.stool$sample.id <- row.names(map.stool)

source_data.1e <- cbind(map.stool$sample.id,  map.stool$DonorID, map.stool$Family.ID, map.stool$Disease.state)

colnames(source_data.1e) <- c("sample.id",  "DonorID", "Family.ID", "Disease.state")

source_data.1e <- cbind("", source_data.1e, "", "", stool.plsda$variates$X)
source_data.1e[1,1] <- "Figure 1e"
source_data.1e[1,7] <- "PLSDA coordinates"

write.csv(source_data.1e,"source_data/source_data_1e.csv", row.names = F)
res <- capture.output(res)
writeLines(res, "Statistical_summaries/Fig_1e_tests.txt")

