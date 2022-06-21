# functions used for Famial Dysauto pub
# Nick Pinkham


######################### FUNCTIONS ##############################3

get.close <- function(map, day = 100){

  map <- map[ !is.na(map$Age.Days) , ]
  mice <- unique(map$Mouse.ID)

  sample.pick <- NULL
  for(i in 1: length(mice)){

    map.i <- map[map$Mouse.ID == mice[i] , ]
    a <- abs(map.i$Age.Days - day)

    sample.pick[i] <- map.i$ID.Gen.Tmt.Age[a == min(a)]
  }

  sample.pick <- sample.pick[!is.na( sample.pick)]

  map.pick <- map[map$ID.Gen.Tmt.Age %in% sample.pick , ]
  map.pick <- map.pick[map.pick$Age.Bin == day , ]
  return(map.pick)
}

norm <- function(samples.pick, metabolites){

  meta.pick <- metabolites[row.names(metabolites) %in% samples.pick , ]

  write.csv(meta.pick, "processing/meta_pick.csv" )

  Metabo_Obj <- InitDataObjects("conc", "stat", FALSE)
  Metabo_Obj <- Read.TextData(Metabo_Obj, "processing/meta_pick.csv")

  #Perform data processing
  Metabo_Obj <- SanityCheckData(Metabo_Obj)
  Metabo_Obj <- ReplaceMin(Metabo_Obj);
  Metabo_Obj <- PreparePrenormData(Metabo_Obj)
  Metabo_Obj <- Normalization(Metabo_Obj, "SumNorm", "LogNorm", "AutoNorm", ref= NULL, ratio = FALSE, ratioNum = 20)
  #Metabo_Obj <- PlotNormSummary(Metabo_Obj, "norm_O_", "png",72, width = NA)
  #Metabo_Obj <- PlotSampleNormSummary(Metabo_Obj, "snorm_O_", "png",72, width = NA)
  Metabo_Obj <- SaveTransformedData(Metabo_Obj)

  meta.normalized <- read.csv("data_normalized.csv", row.names = 1, header = T)

  return(meta.normalized)
}


remove_rare <- function( table , cutoff_pro ) {

  table <- t(table)# transpose to keep "tidy" ; easier that rewriting function...
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) )
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  res <- t(table [ row2keep , , drop=F ])
  res <- as.data.frame(res)

  return(res )
}


pheno.score <- function(map, days ){

  map.pick <- map[!is.na( map$Phenotype.score) & map$Age.Bin == days, ]

  map.pick <- get.close(map.pick, days)

  res <- t.test(map.pick$Phenotype.score ~ map.pick$Treatment.type)

  map.pick$col <- as.factor(map.pick$col)
  vioplot(map.pick$Phenotype.score ~ map.pick$discription,
          main = paste("t test p value =", round(res$p.value, 5),
                       "\n diff between means", round(diff(res$estimate), 5),
                       "\n ~",days, " days"),
          ylab = "Phenotype Score", xlab = "",
          col = as.character(unique(map.pick$col)),
          ylim = c(0,10))
}

pairwiseRF <- function(mat, map, dis_1 = "T 0", dis_2 = "B 7",
                       num_of_col = 50, Accuray.cutoff = 0,
                       Factor = 1, plot = T){
  # report mats that are different between two treatments based on RF analysis
  set.seed(42)
  mat.pick <- mat[map$discription %in% c(dis_1, dis_2 ), ]
  map.pick <- map[map$discription %in% c(dis_1, dis_2 ), ]

  mat.pick$Group <- NULL

  mat.pick.rare_removed <- remove_rare(table = mat.pick, cutoff_pro = 0.2)

  spliter <- as.factor(as.character(map.pick$discription))

  fit <- randomForest::randomForest(x = mat.pick.rare_removed,
                                    y = spliter,
                                    importance=TRUE,
                                    proximities=TRUE,
                                    ntree = 5001)
  RF.sig <- rfUtilities::rf.significance(x = fit,
                                         xdata = mat.pick.rare_removed,
                                         nperm = 101,
                                         ntree = 501)
  print(fit)
  print(RF.sig)

  if(RF.sig$pValue < 0.05){

    fit_imp <- as.data.frame( fit$importance )
    fit_imp$features <- rownames(fit_imp )
    fit_imp_sorted <- arrange( fit_imp  , desc(fit_imp$MeanDecreaseGini)  )
    colors <- vector(length = ncol(mat.pick.rare_removed))

    for(j in 1:ncol(mat.pick.rare_removed)){
      i <- fit_imp_sorted$features[j]
      t1.mean <- mean(mat.pick[map.pick$discription == dis_1, which(colnames(mat.pick) == i)])
      t2.mean <- mean(mat.pick[map.pick$discription == dis_2, which(colnames(mat.pick) == i)])
      if( t1.mean >  t2.mean){
        colors[j] <- "cadetblue3"
      }else{
        colors[j] <- "orange3"
      }
    }

    if(plot){
      barplot(fit_imp_sorted$MeanDecreaseGini[1:num_of_col],
              names.arg = fit_imp_sorted[1:num_of_col,"features"],
              ylab= "Mean Decrease in Accuracy (Variable Importance)",
              las = 2,
              col = colors,
              main= paste(dis_1,"or", dis_2, "Classification RF"))

      legend("topright",
             legend = paste("elivated in", c(dis_1, dis_2)),
             fill = c("cadetblue3", "orange3"))

    }
    res <- fit_imp_sorted$MeanDecreaseGini
    names(res) <- fit_imp_sorted$features

    if(plot){
      abline(h = median(res) + sd(res) * Factor,
             col = 2, lty = 2)
      text(50, median(res) + sd(res) * Factor + .01, paste(" median Feature impotance +", Factor, "sd"),
           col = 2, font = 2)
    }
    res.pick <- res[res > median(res) + (sd(res) * Factor)]
    colors.pick <- colors[res > median(res) + (sd(res) * Factor)]

    mat.bloom <- res.pick[colors.pick == "orange3"]
    mat.bust <- res.pick[colors.pick == "cadetblue3"]

    Result <- list(  mat.bust, mat.bloom)
    names(Result) <- c(dis_1, dis_2)

    return(Result)
  }else{
    return("VARIABLE p val > 0.05")
  }
}

#' Plots PCoA of metabolite  euclidean dissimilarity and colors points  by group
#'
#'
#' @param infile map and metabolite matrix
#' @return PERMANOVA result
plot.meta <- function(meta.pick, map.pick, return.points = F){

  set.seed(42)
  res <- vegan::adonis(dist(meta.pick) ~ map.pick$discription, permutations = 9999)
  print(res)
  res <- res$aov.tab

  ord <- pco(dist(meta.pick))

  plot(ord$points,
       pch = 21,
       col = 8,
       cex = 2,
       bg = map.pick$col,
       ylim = c(min(ord$points[,2]) * 1.1 - .65, max(ord$points[,2]) * 1.1),
       xlab = paste0("PCO 1\n(",
                     round(ord$eig[1] / sum(ord$eig) * 100, 2),
                     "%)"),
       ylab = paste0("PCO 2 (",
                     round(ord$eig[2] / sum(ord$eig) * 100, 2),
                     "%)"))
  vegan::ordispider(ord, map.pick$discription, col = 8)

  legend("bottomright",
         fill = unique(map.pick$col),
         legend = unique(map.pick$discription))

  res1 <- c( round(res$`Pr(>F)`[1], 4),
             round(res$F[1], 3),
             nrow(map.pick))
  names(res1) <- c("p.val", "f.stat", "n")

  title(paste("PERMANOVA\np val =", res1[1],
              "\nf stat =", res1[2],
              "\nn=", res1[3]))
  if(return.points){

    return(ord$points)

  }else{

    return(res1)
  }
}





group.nmds <- function(otu, predictor, groups, col){
  # predictor is a vector of the catagorical predicting variable
  # groups are the paticular catagories you want to explore
  # col is a vertor coresponding to the groups in the predictor

  all.groups <-  unique(predictor)
  others <- all.groups [!all.groups %in% groups]
  others <- paste(  others , collapse = "\nand ")

  nmds <- metaMDS(otu, k = 2, distance = "bray")
  plot(nmds$points, pch = 19, col = 8, cex  = 2,
       xlim = c(min(nmds$points) * 1.1, max(nmds$points) * 1.5))


  nmds.pick <- nmds$points[predictor %in% groups , ]
  col <- col[ predictor %in% groups]
  predictor <- predictor[ predictor %in% groups]

  ordispider(nmds.pick,
             predictor,
             col = unique(col)[order(unique(predictor))],
             label = F, lwd = 2)

  points(nmds.pick, pch = 21, bg = col, cex = 2)


  legend("topright", fill = c(unique(col), 8),
         legend = c(unique(predictor), others))


}



group.nmds2 <- function(otu, predictor, groups, col){
  # predictor is a vector of the catagorical predicting variable
  # groups are the paticular catagories you want to explore
  # col is a vertor coresponding to the groups in the predictor

  otu <- otu[predictor %in% groups , ]
  col <- col[ predictor %in% groups]
  predictor <- predictor[ predictor %in% groups]


  nmds <- metaMDS(otu, k = 2, distance = "bray")
  plot(nmds$points, pch = 19, col = col, cex  = 2)
  ordispider(nmds$points,   predictor,
             col = unique(col)[order(unique(predictor))],
             label = F, lwd = 2)
  legend("topright", fill = unique(col),
         legend = unique(predictor))
}

paired.map <- function(map, pats){

  pat.id <- rel.id <- NULL

  for(i in 1: length(pats)){

    fam <- map$Family.ID[map$DonorID == pats[i]]
    fam.i <- map[map$Family.ID == fam[1] , ]

    pat.i <- map[map$DonorID == pats[i] , ]
    pat.i <- pat.i[ sample(x = 1 : nrow(pat.i), size = 1) , ]
    pat.id[i] <-   pat.i$sample.id

    rel.i <- fam.i[fam.i$Disease.state == "Relative" , ]
    rel.i <- rel.i[abs(pat.i$Collection.date - rel.i$Collection.date) < 90 , ]


    if(nrow(rel.i) > 0){

      rel.id[i] <- rel.i$sample.id[sample(x = 1 : nrow(rel.i), size = 1)]
    }else{
      rel.id[i] <- NA
    }
  }

  map.pat <- map[match(pat.id, map$sample.id) , ]
  map.rel <- map[match(rel.id, map$sample.id) , ]

  map.pick <- rbind(map.pat, map.rel)

  return(  map.pick)

}


pair.meta <- function(map.meta, meta, comps){

#  par(mfrow=c(3,1))

  res <- NULL
  tp <- sort(unique(map.meta$Age.Bin))

  for(i in 1 : length(tp)){


    ################ get the samples we are working with ###################
    samples.pick <- map.meta$ID.Gen.Tmt.Age[map.meta$Age.Bin == tp[i]
                                            & map.meta$discription
                                            %in% comps]

    meta.pick <- norm(samples.pick, meta)
    meta.pick$Label <- NULL

    print("metabolites were normalized")

    map.pick <- map.meta[map.meta$ID.Gen.Tmt.Age %in% samples.pick , ]
    map.pick <- map.pick[ match(row.names(meta.pick),
                                map.pick$ID.Gen.Tmt.Age) , ]


    ###### PERMANOVA ######
    par(mar=c(5.1, 4.1, 7.1, 2.1))

    res1 <- adonis(dist(meta.pick) ~ map.pick$discription, permutations = 9999)

    res1 <- res1$aov.tab

    ord <- pco(dist(meta.pick))

    plot(ord$points,
         pch = 21,
         bg = map.pick$col,
         cex = 2,
         ylim = c(min(ord$points[,2]) * 1.35, max(ord$points) * 1.1),
         xlab = paste0("PCO 1\n(",
                       round(ord$eig[1] / sum(ord$eig) * 100, 2),
                       "%)"),
         ylab = paste0("PCO 2 (",
                       round(ord$eig[2] / sum(ord$eig) * 100, 2),
                       "%)"))
    ordispider(ord, map.pick$discription, col = unique(map.pick$col))

    legend("bottom",
           fill = unique(map.pick$col),
           legend = unique(map.pick$discription))

    res1 <- c( round(res1$`Pr(>F)`[1], 4),
               round(res1$F[1], 3),
               nrow(map.pick))
    names(res1) <- c("p.val", "f.stat", "n")

    title(paste("PERMANOVA\np val =", res1[1],
                "\nf stat =", res1[2],
                "\nn=", res1[3]))
    legend("top", legend = paste("day", tp[i]),  bty="n" )

    res <- rbind(res, res1)
  }


  group.1 <- comps[1]
  group.2 <- comps[2]

  res <- cbind(group.1, group.2, res)

  return(res)
}

pair.meta2 <- function(map.meta, meta, comps){

  l.1 <- NULL
  l.2 <- NULL

  l.1[[1]] <- list()
  l.2[[1]] <- list()

  tp <- sort(unique(map.meta$Age.Bin))

  for(i in 1 : length(tp)){


    ################ get the samples we are working with ###################
    samples.pick <- map.meta$ID.Gen.Tmt.Age[map.meta$Age.Bin == tp[i]
                                            & map.meta$discription
                                            %in% comps]

    meta.pick <- norm(samples.pick, meta)
    meta.pick$Label <- NULL

    print("metabolites were normalized")

    map.pick <- map.meta[map.meta$ID.Gen.Tmt.Age %in% samples.pick , ]
    map.pick <- map.pick[ match(row.names(meta.pick),
                                map.pick$ID.Gen.Tmt.Age) , ]


    ###### PERMANOVA ######
    par(mar=c(5.1, 4.1, 7.1, 2.1))

    set.seed(42)
    res1 <- adonis(dist(meta.pick) ~ map.pick$discription, permutations = 9999)

    res1 <- res1$aov.tab

    ord <- pco(dist(meta.pick))

    plot(ord$points,
         pch = 21,
         bg = map.pick$col,
         cex = 2,
         ylim = c(min(ord$points[,2]) * 1.35, max(ord$points) * 1.1),
         xlab = paste0("PCO 1\n(",
                       round(ord$eig[1] / sum(ord$eig) * 100, 2),
                       "%)"),
         ylab = paste0("PCO 2 (",
                       round(ord$eig[2] / sum(ord$eig) * 100, 2),
                       "%)"))
    ordispider(ord, map.pick$discription, col = unique(map.pick$col))


    legend("bottom",
           fill = unique(map.pick$col),
           legend = unique(map.pick$discription))

    res1 <- c( round(res1$`Pr(>F)`[1], 4),
               round(res1$F[1], 3),
               nrow(map.pick))
    names(res1) <- c("p.val", "f.stat", "n")

    title(paste("PERMANOVA\np val =", res1[1],
                "\nf stat =", res1[2],
                "\nn=", res1[3]))
    legend("top", legend = paste("day", tp[i]),  bty="n" )

    ##### RF #####

    rf.res <- pairwiseRF(mat = meta.pick,
                         map = map.pick,
                         dis_1 = unique(map.pick$discription)[1],
                         dis_2 = unique(map.pick$discription)[2],
                         plot = F)

    if(   rf.res[1] != "VARIABLE p val > 0.05" & res1[1] < 0.05){

      print("RUNNING RANDOM FOREST")

      aa <- match(names(rf.res[[1]]), colnames(meta.pick))
      bb <- match(names(rf.res[[2]]), colnames(meta.pick))

      meta.dif <- c(aa, bb)

      m1 <- meta.pick[map.pick$discription == comps[1],
                      meta.dif]

      m2 <- meta.pick[map.pick$discription == comps[2],
                      meta.dif]

      a <- 0 : (length(m1) -1) * 4

      par(mar=c(4.1, 12.1, 6.1, 4.1))

      boxplot(m1, las = 2, horizontal = TRUE, at = a,
              xlim = c(-.5, max(a)+2),
              col = map.pick$col[map.pick$discription == comps[1]])
      boxplot(m2, las = 2, horizontal = TRUE, yaxt="none",
              at = a + 1, add= T, names = rep("", length(a)),
              col = map.pick$col[map.pick$discription == comps[2]])

      abline(h = length(aa) * 4 - 1.5, lty = 1, col = 3, lwd = 2)

      l.1[[i]] <- names(rf.res[[1]])
      l.2[[i]] <- names(rf.res[[2]])

      meta.pick2 <- meta.pick[,-c(aa, bb)]
      meta.pick2 <- cbind(map.pick$discription, meta.pick2)

      meta.pick2 <- norm(map.pick$ID.Gen.Tmt.Age, meta.pick2)
      meta.pick2$Label <- NULL

      res2 <- adonis(dist(meta.pick2) ~ map.pick$discription, permutations = 9999)
      res2 <- res2$aov.tab
      res2 <- c( round(res2$`Pr(>F)`[1], 4),
                 round(res2$F[1], 3),
                 nrow(map.pick))
      names(res2) <- c("p.val", "f.stat", "n")

      res3 <- res1 - res2

      if(res2[1] < .1){
      title(paste("Metabolites identified as significant by RF",
                  "\nPERMANOVA after metabolites removed\np val =", res2[1],
                  "f stat reduction =", res3[2]))
      }else{
        title(paste("Metabolites identified as significant by RF",
                    "\nPERMANOVA after metabolites removed\np val =", res2[1]))
      }

    }else{
      plot.new()
      text(.5,.5, "Random Forrest\np val > 0.05")
    }
  }

  res <- list(l.1, l.2)
  names(res) <- c(comps[1], comps[2])

  return(res)

}

pair.otu <- function(map.otu, otu, comps, plot.bar = T){

  l.1 <- NULL
  l.2 <- NULL

  l.1[[1]] <- list()
  l.2[[1]] <- list()

  tp <- sort(unique(map.otu$Age.Bin))

  for(i in 1 : length(tp)){


    ################ get the samples we are working with ###################
    map.pick <- map.otu[map.otu$Age.Bin == tp[i] , ]
    #    & map.otu$discription %in% comps , ]
    otu.pick <- otu[map.otu$Age.Bin == tp[i] , ]
    #   & map.otu$discription %in% comps , ]

    map.pick <- get.close(map = map.pick, day = tp[i])
    otu.pick <- otu[ match(map.pick$X16SFile.ID, row.names(otu)), ]

    group.nmds(otu.pick, predictor =  map.pick$discription,
               groups = comps,
               col = map.pick$col)

    map.pick2 <- map.pick[map.pick$discription %in% comps  , ]
    otu.pick2 <- otu.pick[map.pick$discription %in% comps  , ]

    set.seed(42)
    res <- adonis(vegdist(otu.pick2) ~  map.pick2$discription, permutations = 9999)

    res1 <- c( round(res$aov.tab$`Pr(>F)`[1], 4),
               round(res$aov.tab$F[1], 3),
               nrow(map.pick2))
    names(res1) <- c("p.val", "f.stat", "n")

    title(paste(names(res1), "=", res1))

    ##### RF #####

    rf.res <- pairwiseRF(mat = otu.pick2,
                         map = map.pick2,
                         dis_1 = unique(map.pick2$discription)[1],
                         dis_2 = unique(map.pick2$discription)[2],
                         plot = F)


    if(   rf.res[1] != "VARIABLE p val > 0.05"){

      aa <- match(names(rf.res[[1]]), colnames(otu.pick2))
      bb <- match(names(rf.res[[2]]), colnames(otu.pick2))

      meta.dif <- c(aa, bb)

      m1 <- otu.pick2[map.pick2$discription == comps[1],
                      aa]

      m2 <- otu.pick2[map.pick2$discription == comps[2],
                      aa]

      a <- 0 : (length(m1) -1) * 4

      par(mar=c(4.1, 6.1, 6.1, 2.1))


      m1 <- sweep(m1, 2 ,colSums(rbind(m1, m2)),`/`)
      m2 <- sweep(m2, 2 ,colSums(rbind(m1, m2)),`/`)

      if(plot.bar){

        boxplot(m1[,rev(order(colnames(m1)))], las = 2, horizontal = TRUE,
                at = a,
                ylim = range(c(m1, m2)),
                xlab = "Proportion of total counts of OTU\n(within comparison)",
                col = map.pick$col[map.pick$discription == comps[1]])

        boxplot(m2[,rev(order(colnames(m1)))] , las = 2, horizontal = TRUE, yaxt="none",
                at = a + 1, add= T, names = rep("", length(a)),
                col = map.pick$col[map.pick$discription == comps[2]])


        otu.pick2 <- otu.pick[,-c(aa, bb)]

        otu.pick2 <- rrarefy(otu.pick2, min(rowSums(otu.pick2)))

        res2 <- adonis(otu.pick2 ~ map.pick$discription, permutations = 9999)
        res2 <- res2$aov.tab
        res2 <- c( round(res2$`Pr(>F)`[1], 4),
                   round(res2$F[1], 3),
                   nrow(map.pick))
        names(res2) <- c("p.val", "f.stat", "n")

        res2

        if(res2[1] < .1){
          title(paste("Metabolites identified as significant by RF",
                      "\nPERMANOVA after metabolites removed\np val =", res2[1],
                      "f stat =", res2[2]))
        }else{
          title(paste("Metabolites identified as significant by RF",
                      "\nPERMANOVA after metabolites removed\np val =", res2[1]))
        }


        m1 <- otu.pick2[map.pick2$discription == comps[1],
                        bb]

        m2 <- otu.pick2[map.pick2$discription == comps[2],
                        bb]


        m1 <- sweep(m1, 2 ,colSums(rbind(m1, m2)),`/`)
        m2 <- sweep(m2, 2 ,colSums(rbind(m1, m2)),`/`)

        a <- 0 : (length(m1) -1) * 4

        boxplot(m1[,rev(order(colnames(m1)))], las = 2, horizontal = TRUE,
                at = a,
                ylim = range(c(m1, m2)),
                xlab = "Proportion of total counts of OTU\n(within comparison)",
                col = map.pick$col[map.pick$discription == comps[1]])

        boxplot(m2[,rev(order(colnames(m1)))] , las = 2, horizontal = TRUE, yaxt="none",
                at = a + 1, add= T, names = rep("", length(a)),
                col = map.pick$col[map.pick$discription == comps[2]])

        #abline(h = length(aa) * 4 - 1.5, lty = 1, col = 3, lwd = 2)

        title(paste(comps[2], "\nOTUs identified as significant\nthrough RF analysis"))
      }
      l.1[[i]] <- names(rf.res[[1]])
      l.2[[i]] <- names(rf.res[[2]])

    }else{
      plot.new()
    }
  }
  res <- list(l.1, l.2)
  names(res) <- c(comps[1], comps[2])

  return(res)
}

pair.otu2 <- function(map.otu, otu, comps, plot.bar = T){

  tp <- sort(unique(map.otu$Age.Bin))
  res <- NULL

  for(i in 1 : length(tp)){


    ################ get the samples we are working with ###################
    map.pick <- map.otu[map.otu$Age.Bin == tp[i] , ]
    #    & map.otu$discription %in% comps , ]
    otu.pick <- otu[map.otu$Age.Bin == tp[i] , ]
    #   & map.otu$discription %in% comps , ]

    map.pick <- get.close(map = map.pick, day = tp[i])
    otu.pick <- otu[ match(map.pick$X16SFile.ID, row.names(otu)), ]

    group.nmds(otu.pick, predictor =  map.pick$discription,
               groups = comps,
               col = map.pick$col)

    map.pick2 <- map.pick[map.pick$discription %in% comps  , ]
    otu.pick2 <- otu.pick[map.pick$discription %in% comps  , ]

    set.seed(42)
    res.i <- adonis(vegdist(otu.pick2) ~  map.pick2$discription, permutations = 9999)

    res.i <- c(tp[i],
               round(res.i $ aov.tab$`Pr(>F)`[1], 4),
               round(res.i $ aov.tab$F[1], 3),
               nrow(map.pick2))

    names(res.i) <- c("day", "p.val", "f.stat", "n")

    title(paste(names(res.i), "=", res.i))

    res <- rbind(res, res.i)
  }

  group.1 <- comps[1]
  group.2 <- comps[2]

  res <- cbind(group.1, group.2, res)

  return(res)
}

test.time.otu <- function(group){

  map.pick <- map.otu[map.otu$discription == group, ]
  map.pick <- map.pick[order(map.pick$Age.Bin) , ]
  pick1 <- get.close(map.pick, 100)$X16SFile.ID
  pick2 <- get.close(map.pick, 200)$X16SFile.ID
  pick3 <- get.close(map.pick, 300)$X16SFile.ID

  map.pick2 <- map.otu[map.otu$X16SFile.ID %in% c(pick1, pick3) , ]
  otu.pick2 <- otu[match(map.pick2$X16SFile.ID, rownames(otu)) , ]

  map.pick3 <- map.otu[map.otu$X16SFile.ID %in% c(pick1,
                                                  pick2,
                                                  pick3) , ]
  otu.pick3 <- otu[match(map.pick3$X16SFile.ID, rownames(otu)) , ]

  res <- adonis(vegdist(otu.pick2) ~ map.pick2$Age.Bin)

  print(res)


  col.order <-
    unique(map.pick3$Age.Col)[order(unique(as.character(map.pick3$Age.Bin)))]

  boxplot(  map.pick3$invsimp ~ map.pick3$Age.Bin, col = col.order)

  nmds <- metaMDS(otu.pick3, k = 2, distance = "bray")

  plot(nmds$points, pch = 19,
       col = map.pick3$Age.Col,
       cex  = 2)

#  text(nmds$points[,1], nmds$points[,2], map.pick3$Age.Days, pos = 2)

  ordispider(nmds,
             map.pick3$Age.Bin,
             label = F,
             col = col.order)

  res <- c( round(res$aov.tab$`Pr(>F)`[1], 4),
            round(res$aov.tab$F[1], 3),
            nrow(map.pick2))
  names(res) <- c("p.val", "f.stat", "n")

  title(paste(group, "first age bin to last \nPERMANOVA p val =", res[1],
              "f stat =", res[2],
              "n=", res[3]))

  legend("bottomleft",
         fill = unique(map.pick3$Age.Col),
         legend = unique(map.pick3$Age.Bin))
}

test.time.meta <- function(group){

  map.pick <- map.meta[map.meta$discription == group, ]
  pick1 <- get.close(map.pick, 100)$ID.Gen.Tmt.Age
  pick2 <- get.close(map.pick, 200)$ID.Gen.Tmt.Age
  pick3 <- get.close(map.pick, 300)$ID.Gen.Tmt.Age

  map.pick2 <- map.meta[map.meta$ID.Gen.Tmt.Age %in% c(pick1, pick3) , ]

  meta$Group <- paste0(meta$Group, map.meta$Age.Bin)

  meta.pick <- norm(samples.pick = map.pick2$ID.Gen.Tmt.Age, metabolites = meta)

  meta.pick$Label <- NULL
  res <- adonis(dist(meta.pick) ~ map.pick2$Age.Bin)

  print(res)

  ord <- pco(dist(meta.pick))

  col.order <-
    unique(map.pick2$Age.Col)[order(unique(as.character(map.pick2$Age.Bin)))]

  plot(ord$points,
       pch = 21,
       bg = map.pick2$Age.Col,
       cex = 2,
       ylim = c(min(ord$points[,2]) * 1.35, max(ord$points) * 1.1),
       xlab = paste0("PCO 1\n(",
                     round(ord$eig[1] / sum(ord$eig) * 100, 2),
                     "%)"),
       ylab = paste0("PCO 2 (",
                     round(ord$eig[2] / sum(ord$eig) * 100, 2),
                     "%)"))

  ordispider(ord$points,
             map.pick2$Age.Bin,
             col =   col.order)

  legend("bottom",
         fill =   col.order,
         legend = unique(map.pick2$Age.Bin))

  res <- c( round(res$aov.tab$`Pr(>F)`[1], 4),
             round(res$aov.tab$F[1], 3),
             nrow(map.pick2))
  names(res) <- c("p.val", "f.stat", "n")

  title(paste(group, "First to Last age bin\nPERMANOVA p val =", res[1],
              "f stat =", res[2],
              "n=", res[3]))
}

group.nmds2 <- function(otu, predictor, groups, col){
  # predictor is a vector of the catagorical predicting variable
  # groups are the paticular catagories you want to explore
  # col is a vertor coresponding to the groups in the predictor

  all.groups <-  unique(predictor)
  others <- all.groups [!all.groups %in% groups]
  others <- paste(  others , collapse = "\nand ")

  nmds <- metaMDS(otu, k = 2, distance = "bray")
  plot(nmds$points, pch = 19, col = 8, cex  = 2,
       xlim = c(min(nmds$points) * 1.1, max(nmds$points) * 1.5))


  nmds.pick <- nmds$points[predictor %in% groups , ]
  col <- col[ predictor %in% groups]
  predictor <- predictor[ predictor %in% groups]

  col.map <- aggregate(map.pick$col,
                       by = list(map.pick$discription),
                       FUN = "unique")

  ordispider(nmds.pick,
             predictor,
             col = col.map$x,
             #col = unique(col)[order(unique(predictor))],
             label = F, lwd = 2)

  points(nmds.pick, pch = 21, bg = col, cex = 2)


  legend("topright", fill = c(col.map$x, 8),
         legend = c(col.map$Group.1, others))
}


group.medion <- function(otu, map, group, day){

  map.p <- get.close(map = map.otu[map$discription == group , ], day)
  otu.p <- otu[ match(map.p$X16SFile.ID, row.names(otu)), ]
  med <- apply(otu.p, 2, median)

  return(med)
}

beta.dis <- function(otu, map){

  beta <- betadisper(vegdist(otu), map$Disease.Bin)

  set.seed(999)

  res1 <- adonis(otu ~ map$Disease.Bin)

  res1 <- res1$aov.tab

  res1 <- c( round(res1$`Pr(>F)`[1], 4),
             round(res1$F[1], 3),
             nrow(map))
  names(res1) <- c("p.val", "f.stat", "n")

  plot(beta,
       main = paste("PERMANOVA\np val =", res1[1],
                    "\nf stat =", res1[2],
                    "\nn=", res1[3]))

  text(beta$vectors[,1:2], labels = map$Phenotype.score, pos = 2)


  boxplot(beta$distances ~ map$Disease.Bin)
  t.res <- t.test(beta$distances ~ map$Disease.Bin)
  title(paste("T test p value =", round(t.res$p.value, 3)))

  plot(map$Phenotype.score, beta$distances,
       xlab = "Phenotype score",
       ylab = "BC to distances centroid",
       pch = 21, bg = 8)
  abline(lm( beta$distances ~ map$Phenotype.score), lty = 2, col = 2)
  c.res <- cor.test(map$Phenotype.score, beta$distances)

  title(paste("Pearson's correlation\np val =", round(c.res$p.value, 3),
              "\ncor =", round(c.res$estimate, 3)))

  return(c.res)
}


d2med <- function(otu, map, mediod) {

  otu.comp <- rbind(mediod, otu)
  dist <- as.matrix(vegdist(otu.comp))
  dist2control <- dist[ -1, 1]
  plot(map$Phenotype.score, dist2control,
       xlab = "Phenotype score",
       ylab = "BC to distances centroid",
       pch = 21, bg = 8)
  abline(lm( dist2control ~ map$Phenotype.score), lty = 2, col = 2)

  res <- cor.test(map$Phenotype.score, dist2control)

  title(paste("p val =", res$p.value,
              "\ncor =", res$estimate))
}


