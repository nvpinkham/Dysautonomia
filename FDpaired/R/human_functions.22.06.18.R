# functions used for Famial Dysauto pub
# Nick Pinkham 

library(vegan)
library(labdsv)
library(cluster)
library(vioplot)
library(gplots)
library(RColorBrewer) 
library(dendextend)
library(dplyr)
library(mixOmics)
library(vegan3d)
library(scatterplot3d)
library(lme4)
options(scipen=999)


rs <-list.files("/Users/nickpinkham/Desktop/Nature2022_submission/data/MetaboAnalystR-master/R", pattern = ".R", full.names = T)

for(i in rs[-1]){
  # Load MetaboAnalystR functions
    source(i)
}


dex <- function(x){
  x <- gsub("[.]", "-", x)
  b <- grep("X.-", x)
  c <- x[b] 
  x[b] <- substr(c, 2, 100)
  return(x)
}

group.nmds <- function(otu, predictor, groups, col){
  # predictor is a vector of the catagorical predicting variable
  # groups are the paticular catagories you want to explore
  # col is a vertor coresponding to the groups in the predictor
  
  all.groups <-  unique(predictor)
  others <- all.groups [!all.groups %in% groups]
  others <- paste(  others , collapse = "\nand ")
  
  nmds <- metaMDS(otu, k = 2, distance = "bray")
  plot(nmds$points, pch = 19, col = 8, cex  = 2)
  
  
  nmds.pick <- nmds$points[predictor %in% groups , ]
  col <- col[ predictor %in% groups]
  predictor <- predictor[ predictor %in% groups]
  
  ordispider(nmds.pick,   predictor, 
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

#' Generate a map with only paired samples 
#'
#' Samples are paired randomly sith sample function. 
#'
#' @param infile original map and list of patients
#' @return map with random pairs
#' @export

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


#' Gets the sample closest to a day 
#'
#'
#' @param infile map and day as integer
#' @return 
get.close <- function(map, day){
  
  mice <- unique(map$Mouse.ID)
  mice
  
  sample.pick <- NULL
  for(i in 1: length(mice)){
    
    map.i <- map[map$Mouse.ID == mice[i] , ]
    a <- abs(map.i$age - day)
    
    sample.pick[i] <- map.i$Sample.ID[a == min(a)]
  }
  
  map.pick <- map[map$Sample.ID %in% sample.pick , ]
  return(map.pick)
}

#' Normalizes metabolite matrix with metaboanalyst
#'
#'
#' @param infile list of samples being used and matrix of raw metabolite concentrations 
#' @return glog normalized metabolite concentrations 
#' 
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
  
  file.remove("data_original.csv")
  file.remove("data_normalized.csv")
  file.remove("data_processed.csv")
  sapply(list.files(pattern = ".qs"), file.remove)
  
  return(meta.normalized)
}

colMedians <- function(mat){
  res <- apply(mat, 2, median)
  return(res)
}


norm.human <- function(samples.pick, metabolites){
  
  meta.pick <- metabolites[row.names(metabolites) %in% samples.pick , ]
  
  if(!dir.exists("processing")){
    dir.create("processing")
  }
  
  write.csv(meta.pick, "processing/meta_pick.csv" )
  
  Metabo_Obj <- InitDataObjects("conc", "stat", FALSE)
  Metabo_Obj <- Read.TextData(Metabo_Obj, "processing/meta_pick.csv")
  
  #Perform data processing
  Metabo_Obj <- SanityCheckData(Metabo_Obj)
  Metabo_Obj <- ReplaceMin(Metabo_Obj);
  Metabo_Obj <- PreparePrenormData(Metabo_Obj)
  Metabo_Obj <- Normalization(Metabo_Obj, "NULL", "LogNorm", "AutoNorm", ref= NULL, ratio = FALSE, ratioNum = 20)
  #Metabo_Obj <- PlotNormSummary(Metabo_Obj, "norm_O_", "png",72, width = NA)
  #Metabo_Obj <- PlotSampleNormSummary(Metabo_Obj, "snorm_O_", "png",72, width = NA)
  Metabo_Obj <- SaveTransformedData(Metabo_Obj)
  
  meta.normalized <- read.csv("data_normalized.csv", row.names = 1, header = T)
  
  file.remove("data_original.csv")
  file.remove("data_normalized.csv")
  file.remove("data_processed.csv")
  sapply(list.files(pattern = ".qs"), file.remove)
  
  return(meta.normalized)
}


colMedians <- function (x){
  
  res <- apply(x, 2, median)
  return(res)
}


plot.pval <- function(mat, cols, name,  type = "stool",  b4 = "after FDR correction", print = T){
  par(mar=c(10.1,4.1,4.1,2.1))
  par(cex.axis=0.5)
  
  aa <- colMedians(mat)
  
  colfunc <- colorRampPalette(cols)
  
  bb <- colfunc(100)
  cc <- 1 : 100 / 100
  ff <- NULL
  
  for(i in 1 : length(aa)){
    
    dd <- abs(aa[i] - cc)
    ff[i] <- bb[which(dd ==  min(dd))]
  }
  
  if(print){
    setEPS()       # Set postscript arguments
    postscript(paste0(name,
                      Sys.Date(), ".eps"),
               width = 20, height = 10) 
  }
  layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
  par(mar=c(12.1,4.1,4.1,2.1))
  
  
  boxplot(mat, las = 2, ylab = paste("p value", b4), 
          cex = .5,
          pch = 21,
          bg = "lightgrey",
          main = paste(type, "METABOLITES 999 PERMUTATIONS"), col = ff)
  
  abline(h = 0.05, col = 2, lty = 2)
  
  legend_image <- as.raster(matrix(colfunc(100), ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'mean p value', )
  text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
  rasterImage(rev(legend_image), 0, 0, 1,1)
  rect(0,1,1,0)
  
  segments(0,0,1,0, col = 2, lty = 2)
  if(print){
    
    dev.off()
  }
}

plot.effect <- function(mat, cols, name, type = "stool", print = T){
  
  par(mar=c(10.1,4.1,4.1,2.1))
  par(cex.axis=0.5)
  
  aa <- colMedians(mat)
  
  colfunc <- colorRampPalette(cols)
  
  bb <- colfunc(201)
  
  cc <- 0 : 200 / 100 - 1
  ff <- NULL
  
  for(i in 1 : length(aa)){
    
    dd <- abs(aa[i] - cc)
    ff[i] <- bb[which(dd ==  min(dd))]
  }
  
  if(print){
    setEPS()       # Set postscript arguments
    postscript(paste0(name,
                      Sys.Date(), ".eps"),
               width = 20, height = 10) 
  }
  
  layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
  par(mar=c(12.1,4.1,4.1,2.1))
  
  boxplot(mat, las = 2, 
          cex = .5,
          pch = 21,
          ylim = c(-1.1, 1.1),
          ylab = "effect size", 
          bg = "lightgrey",
          main = paste(type, "METABOLITES 999 PERMUTATIONS"), col = ff)
  
  abline(h = 0, col = 2, lty = 2)
  
  legend_image <- as.raster(matrix(colfunc(100), ncol=1))
  plot(c(0,2),c(-1,1),type = 'n', axes = F,xlab = '', ylab = '',
       main = 'effect size', )
  
  text(x=1.5, y = seq(-1,1,l=5), labels = seq(-1,1,l=5))
  rasterImage(rev(legend_image), 0, -1, 1,1)
  rect(0,-1,1,1)
  
  segments(0,0,1,0, col = 2, lty = 2)
  
  if(print){
    
    dev.off()
  }
}

tax.shared <- function(map, otu, tax, binary = F){
  
  tax.familes <- sort(unique(tax$family))
  fam.count <- aggregate(t(otu), list(tax$family), sum)
  row.names(fam.count) <- fam.count$Group.1
  fam.count$Group.1 <- NULL
  fam.count <- as.data.frame(t(fam.count))
  
  filter <- NULL
  for(i in 1:length(map$sample.id)){
    # remove samples that can't be paired 
    
    map.i <- map[i,]
    fam.i <- map[map$Family.ID == map.i$Family.ID,]
    
    pairs.pos <- fam.i[fam.i$Disease.state != map.i$Disease.state , ] #possible pairs 
    
    if(nrow(pairs.pos) > 0){
      
      closest.date <- min(abs(pairs.pos$Collection.date - map.i$Collection.date))
      filter[i] <- closest.date < 90
    }else{
      filter[i] <- FALSE
      
    }
  }
  
  map <- map[filter,]
  pats <- unique(map$DonorID[map$Disease.state == "Patient"])
  map.pick <- paired.map(map, pats)
  
  otu.pick <- otu[match(map.pick$sample.id, row.names(otu)) , ]
  fam.pick <- fam.count[match(map.pick$sample.id, row.names(fam.count)) , ]
  
  #### Run the thing #####
  
  pat.unique <- NULL
  rel.unique <- NULL
  
  pat.shared <- NULL
  rel.shared <- NULL
  
  for(j in 1 : length(pats)){
    
    fam.j <- map.pick$Family.ID[map.pick$DonorID == pats[j]]
    
    p <- otu.pick[map.pick$Disease.state == "Patient" & map.pick$Family.ID == fam.j, ]
    r <- otu.pick[map.pick$Disease.state != "Patient" & map.pick$Family.ID == fam.j, ]
    
    if(binary){
      p[p > 0] <- 1
      r[r > 0] <- 1
    }
    
    p.dummy2 <- p.dummy <- p
    p.dummy[colSums(r) > 0] <- 0# shared OTUs removed
    p.dummy2[colSums(  p.dummy) > 0] <- 0# unique OTUs removed
    
    r.dummy2 <-  r.dummy <- r
    r.dummy[colSums(p) > 0] <- 0# shared OTUs removed
    r.dummy2[colSums(  r.dummy) > 0] <- 0# unique OTUs removed

    p.j <- aggregate(t(p.dummy), list(tax$family), sum)
    r.j <- aggregate(t(r.dummy), list(tax$family), sum)
    
    p1.j <- aggregate(t(p.dummy2), list(tax$family), sum)
    r1.j <- aggregate(t(r.dummy2), list(tax$family), sum)
    
    pat.unique <- cbind(pat.unique, p.j[,2])
    rel.unique <- cbind(rel.unique, r.j[,2])
    
    pat.shared <- cbind(pat.shared, p1.j[,2])
    rel.shared <- cbind(rel.shared, r1.j[,2])
  }
  
  tax.res <- cbind(  pat.unique,
                     pat.shared,
                     rel.shared,
                     rel.unique)
  
  rownames(tax.res) <- tax.familes
  colnames(tax.res) <- c(rep("A Pat unique", ncol(pat.unique)), 
                         rep("B Pat shared", ncol(pat.shared)),
                         rep("C Rel shared", ncol(rel.shared)),
                         rep("D Rel unique", ncol(rel.unique)))
  
  tax.res <- tax.res[ order(colSums(fam.count)) , ]
  return(tax.res)
}


filter.pairs <- function(map){
  # remove samples that can't be paired 
  
  filter <- NULL
  for(i in 1 : nrow(map)){
    
    map.i <- map[i,]
    fam.i <- map[map$Family.ID == map.i$Family.ID,]
    
    pairs.pos <- fam.i[fam.i$Disease.state != map.i$Disease.state , ] #possible pairs 
    
    if(nrow(pairs.pos) > 0){
      
      closest.date <- min(abs(pairs.pos$Collection.date - map.i$Collection.date))
      filter[i] <- closest.date < 90
    }else{
      filter[i] <- FALSE
      
    }
  }
  return(filter)
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

#' make violin plot and connects pairs
#'
#' 
#' @return the number of pairs 

paired.violin.beta <- function(map, var, plot.each= F){
  
  set.seed(42)
  map$dummy <- 1
  map$dummy[map$Disease.state != "Patient" ] <- 2
  map$dummy  <- map$dummy + sample(-100 : 100, nrow(map), replace = T)/2000
  
  map.pat <- map[map$Disease.state == "Patient", ]
  map.rel <- map[map$Disease.state == "Relative", ]
  
  var.pat <- var[map$Disease.state == "Patient" ]
  var.rel <- var[map$Disease.state == "Relative"]
  
  num.pairs <- NULL
  
  vioplot(var.pat,  var.rel, drawRect = F, col=unique(map$col), 
          names = c("Patients", "Relatives"))
  
  segments(.7, mean(var.pat),
           1.3, mean(var.pat), lwd = 2)
  
  segments(1.7, mean(var.rel),
           2.3, mean(var.rel), lwd = 2)  
  
  families <- sort(unique(map$Family.ID))

  for(i in 1: length(families)){
  
    map.pat.i <- map.pat[map.pat$Family.ID ==  families[i] , ]
    map.rel.i <- map.rel[map.rel$Family.ID ==  families[i] , ]
    
    var.pat.i <- var.pat[map.pat$Family.ID ==  families[i] ]
    var.rel.i <- var.rel[map.rel$Family.ID ==  families[i] ]
    
    for(j in 1 : nrow(map.pat.i)){
      
      t.diff <- abs(map.pat.i$Collection.date[j] -  map.rel.i$Collection.date) 
      
      if(min(t.diff) <  90){
        
        map.rel.i.j <- map.rel.i[t.diff < 90 , ]
        var.rel.i.j <- var.rel.i[t.diff < 90  ]
        
        
        for(k in 1 :  nrow( map.rel.i.j)){
          
          segments(y0 = var.pat.i[j],  x0 =  map.pat.i$dummy[j], 
                   y1 =  var.rel.i.j[k], x1 =  map.rel.i.j$dummy[k], 
                   col = 8)
          
          num.pairs <- c(num.pairs, 1)
        }
      }
    }
    
    points(map.pat.i$dummy, var.pat.i, pch = 21, bg = 8)
    points(map.rel.i$dummy, var.rel.i, pch = 21, bg = 8)
    
  }
  
  return(  length(num.pairs ))
}


heat <- function(matrix, colors, x.adj = -3){
  
  colfunc <- colorRampPalette(colors)
  
  bb <- colfunc(100)
  cc <- 1 : 100 / 100
  ff <- NULL
  
  matrix.col <- trans01(matrix)
  
  for(i in 1 : nrow(matrix)){
    
    aa <- as.numeric(matrix.col[i,])
    
    column <- row.names(matrix)[i]
    
    text((x.adj + i - .5),
         cex = .75, font = 3, 
         ncol(matrix) * 1.2, 
         paste( paste(rep("  ", nchar(column)), collapse = ""), 
                column), srt = 45)
    
    points((x.adj + i - .5),
           ncol(matrix) * 1.2, 
           pch = "|", bg = 2)
    
    for(j in 1 : ncol(matrix)){
      
      dd <- abs(cc - aa[j])
      col.here <- bb[which(dd ==  min(dd))]
      matrix.col[i, j] <- col.here
      
      rect(xleft = i + x.adj, ybottom = (j * 1.2) - 1.1, 
           xright = (i - 1) + x.adj, ytop = (j * 1.2) + .1, 
           col = col.here)
    }
  }
  return(colfunc)
}

trans01 <- function(x){
  (x-min(x))/(max(x)-min(x))
}



count.pairs <- function(map, var){
  
  map <- map[!is.na(map$Collection.date) , ]
  var <- var[!is.na(map$Collection.date) ]
  
  map.pat <- map[map$Disease.state == "Patient", ]
  map.rel <- map[map$Disease.state == "Relative", ]
  
  var.pat <- var[map$Disease.state == "Patient" ]
  
  
  map.pat <- map.pat[!is.na(var.pat) , ]
  
  families <- sort(unique(map.pat$Family.ID))
  
  num.pairs <- NULL
  
  for(i in 1: length(families)){
    
    map.pat.i <- map.pat[map.pat$Family.ID ==  families[i] , ]
    map.rel.i <- map.rel[map.rel$Family.ID ==  families[i] , ]
    
    for(j in 1 : nrow(map.pat.i)){
      
      t.diff <- abs(map.pat.i$Collection.date[j] -  map.rel.i$Collection.date) 
      map.rel.i.j <- map.rel.i[t.diff < 90 , ]
      
      num.pairs <- c(num.pairs, nrow(  map.rel.i.j ))
      
    }
  }
  return(  sum(num.pairs ))
}


unique.collapse <- function(x){
  paste(sort(unique(x)), collapse = ", ")
}

get.paired <- function(map){
  
  samples <- as.data.frame(map$sample.id)
  samples$`map$sample.id`
  samples$pairs <- NA
  
  for(i in 1: nrow(map)){
    
    map.i.r <- map[ map$Family.ID == map$Family.ID[i] & map$Disease.state == "Relative" , ]
    map.i.r <- map.i.r[ map.i.r$Sample.type == map$Sample.type[i], ]
    
    t.diff <- abs(map$Collection.date[i] - map.i.r$Collection.date) 
    map.i.r.pick <-   map.i.r[which(t.diff < 90), ]
    map.i.r.pick <-  map.i.r.pick[order( map.i.r.pick$sample.id) , ]
    
    samples$paired.samples[i] <- paste(map.i.r.pick$sample.id, collapse = ", ")
    samples$paired.relative[i] <- paste(map.i.r.pick$DonorID, collapse = ", ")
    
    if( map$Disease.state[i] != "Relative" & length(map.i.r.pick$sample.id) > 0 ){
      
      a <- paste(map$sample.id[i],"&", map.i.r.pick$sample.id)
      samples$pairs[i] <- paste(a, collapse = ", ")
      
    }
  }
  return(samples)
}


get.used <- function(map, pick = "Stool 16S"){
  
  pairs.list <- strsplit(map$pairs, ",")
  pairs <- unique(unlist(pairs.list))
  pairs <- pairs[pairs != "NP"]
  
  res <- NULL
  for(i in 1:length(pairs)){
    
    aa <- map$Sample.type[grep(pairs[i], pairs.list)]
    res[i] <- grepl(pick , aa)
  }
  
  res[res] <- "yes"
  res[res == "FALSE"] <- "no"
  
  return(res)
}

