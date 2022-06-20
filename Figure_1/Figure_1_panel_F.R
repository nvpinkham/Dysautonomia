# Figure_1_panel_f

source("human_functions.22.04.11.R")

num.perm = 9

map.stool <- map.stool[!is.na(map.stool$Collection.date) , ]

rs <-list.files("MetaboAnalystR-master/R", pattern = ".R", full.names = T)

for(i in rs[-1]){
  source(i)
}

filter <- NULL
for(i in 1 : nrow(map.stool)){
  # remove samples that can't be paired 
  
  map.i <- map.stool[i,]
  fam.i <- map.stool[map.stool$Family.ID == map.i$Family.ID,]
  
  pairs.pos <- fam.i[fam.i$Disease.state != map.i$Disease.state , ] #possible pairs 
  
  if(nrow(pairs.pos) > 0){
    
    closest.date <- min(abs(pairs.pos$Collection.date - map.i$Collection.date))
    filter[i] <- closest.date < 90
  }else{
    filter[i] <- FALSE
    
  }
}

map.stool  <- map.stool[filter,]
meta.stool <- meta.stool[match(map.stool$sample.id, row.names(meta.stool)) , ]
all(map.stool$sample.id == row.names(meta.stool))
pats <-  unique(map.stool$DonorID[map.stool$Disease.state == "Patient"])

#meta.stool$Label <- NULL
meta.stool$Ethanol <- NULL # stored in EtOH
meta.stool$X <- NULL

stool.t.perm <- matrix(ncol = ncol(meta.stool) -1, nrow = num.perm)
colnames(stool.t.perm) <- colnames(meta.stool[-1])
stool.t.e.perm <- stool.t.perm.no <- stool.t.perm
#eff <- array(dim = c(ncol(meta.stool), length(pats), num.perm))
fold <-  matrix(ncol = ncol(meta.stool) - 1, nrow = num.perm * length(pats))
dim(fold)

d <- length(pats)
permutation <- NULL
comps <- NULL
comps.list <- list()

df <- NULL

map.stool$Collection.date <- as.Date(map.stool$Collection.date, "%m/%d/%y")

for(j in  1 : num.perm){
  print(j)
  
  map.i <- paired.map(map = map.stool, pats)
  
  ## DON'T RESAMPLE Individuals 
  duple.samps <- all(!duplicated(map.i$sample.id))
  while (  !duple.samps ) {
    map.i <- paired.map(map = map.stool, pats)
    duple.samps <- !all(duplicated(map.i$sample.id))
    print("resampleing because one or more samples were used 
          more than once")
  }
  
  p <- match(map.i$sample.id, row.names(meta.stool))
  
  meta.i <- meta.stool[p,]
  map.i.ds <- map.stool$sample.id[p]# if a family has 2 patient they may be paired to the same relative
  
  samps.pick <- map.stool$sample.id[p]
  
  meta.stool.norm <- norm.human(samples.pick = samps.pick, metabolites = meta.stool)
  meta.stool.norm$Label <- NULL
  
  meta.i <- meta.stool.norm[match(map.i$sample.id, row.names(meta.stool.norm)) , ]
  
  for(i in 1 : ncol(meta.stool.norm)){
    
    s <-  meta.i[,i]
    
    res <- t.test(s[map.i$Disease.state == "Patient"],
                  s[map.i$Disease.state != "Patient"], paired = T)
    
    stool.t.perm[j,i] <- res$p.value
    
    #  difs <- meta.stool.norm$Choline[map.i$Disease.state != "Patient"] - meta.stool.norm$Choline[map.i$Disease.state == "Patient"]
    #  stool.t.e.perm[j,i] <- mean(difs)
    
    stool.t.e.perm[j,i] <- res$estimate
    
    here <-  (((j - 1) * d) + 1) : (j * d) 
    
    fold[ here , i] <- s[map.i$Disease.state == "Patient"] /s[map.i$Disease.state != "Patient"]
    
    permutation[here] <- rep(j, d)
    
    comp <- paste(map.i$sample.id[map.i$Disease.state == "Patient"],
                  "to",
                  map.i$sample.id[map.i$Disease.state != "Patient"])
    
    
    comps[here] <- comp
    
    comps.list[[j]] <- comp
  }
  df[j] <- length(s[map.i$Disease.state == "Patient"])
}

stool.fold <- cbind(comps, permutation, fold)

stool.t.e.perm <- stool.t.e.perm[, order(colMedians(stool.t.perm))  ]
stool.t.perm <- stool.t.perm[, order(colMedians(stool.t.perm))  ]

stool.t.fdr <- stool.t.perm

for(i in 1 : nrow(stool.t.perm)){
  
  stool.t.fdr[i,] <- p.adjust(stool.t.perm[i,], method = "fdr")
}

sort(colMeans(stool.t.fdr))

dir.create("Permutation_Result")

write.csv(stool.t.perm, "Permutation_Result/stool_perm_p_vals.csv")
write.csv(stool.t.fdr, "Permutation_Result/stool_perm_p_vals_fdr.csv")

write.csv(stool.t.e.perm, "Permutation_Result/stool_perm_effect_sizes.csv")
write.csv(stool.fold, "Permutation_Result/stool_perm_fold.csv")
################################################################################

filter <- NULL
for(i in 1 : nrow(map.serum)){
  # remove samples that can't be paired 
  
  map.i <- map.serum[i,]
  fam.i <- map.serum[map.serum$Family.ID == map.i$Family.ID,]
  
  pairs.pos <- fam.i[fam.i$Disease.state != map.i$Disease.state , ] #possible pairs 
  
  if(nrow(pairs.pos) > 0){
    
    closest.date <- min(abs(pairs.pos$Collection.date - map.i$Collection.date))
    filter[i] <- closest.date < 90
  }else{
    filter[i] <- FALSE
    
  }
}

map.serum  <- map.serum[filter,]
meta.serum <- meta.serum[match(map.serum$sample.id, row.names(meta.serum)) , ]
all(map.serum$sample.id == row.names(meta.serum))
pats <-  unique(map.serum$DonorID[map.serum$Disease.state == "Patient"])

#meta.serum$Label <- NULL
meta.serum$Ethanol <- NULL # stored in EtOH
meta.serum$X <- NULL

serum.t.perm <-matrix(ncol = ncol(meta.serum) -1, nrow = num.perm)
colnames(serum.t.perm) <- colnames(meta.serum[-1])
serum.t.e.perm <- serum.t.perm.no <- serum.t.perm
#eff <- array(dim = c(ncol(meta.serum), length(pats), num.perm))
fold <-  matrix(ncol = ncol(meta.serum) - 1, nrow = num.perm * length(pats))
dim(fold)

d <- length(pats)
permutation <- NULL
comps <- NULL
comps.list <- list()

df <- NULL

map.serum$Collection.date <- as.Date(map.serum$Collection.date, "%m/%d/%y")

for(j in  1 : num.perm){
  print(j)
  
  map.i <- paired.map(map = map.serum, pats)
  
  ## DON'T RESAMPLE Individuals 
  duple.samps <- all(!duplicated(map.i$sample.id))
  while (  !duple.samps ) {
    map.i <- paired.map(map = map.serum, pats)
    duple.samps <- !all(duplicated(map.i$sample.id))
    print("resampleing because one or more samples were used 
          more than once")
  }
  
  
  
  p <- match(map.i$sample.id, row.names(meta.serum))
  
  meta.i <- meta.serum[p,]
  map.i.ds <- map.serum$sample.id[p]# if a family has 2 patient they may be paired to the same relative
  
  samps.pick <- map.serum$sample.id[p]
  
  meta.serum$Label <- sample(meta.serum$Label)
  meta.serum.norm <- norm.human(samples.pick = samps.pick, metabolites = meta.serum)
  meta.serum.norm$Label <- NULL

  meta.i <- meta.serum.norm[match(map.i$sample.id, row.names(meta.serum.norm)) , ]
  
  for(i in 1 : ncol(meta.serum.norm)){
    
    s <-  meta.i[,i]
    
    res <- t.test(s[map.i$Disease.state == "Patient"],
                  s[map.i$Disease.state != "Patient"], paired = T)
    
    
    meta.raw.i <- meta.serum.norm[match(map.i$sample.id,
                                        row.names(meta.serum.norm)), ]
    
    row.names(meta.raw.i) == map.i$sample.id
    s.raw <- meta.raw.i[,i]
    

    serum.t.perm[j,i] <- res$p.value
    serum.t.e.perm[j,i] <- mean(s.raw[map.i$Disease.state == "Patient"] / s.raw[map.i$Disease.state != "Patient"])
  }
  df[j] <- length(s[map.i$Disease.state == "Patient"])
  
}

fold <- cbind(comps, permutation, fold)

serum.t.e.perm <- serum.t.e.perm[, order(colMedians(serum.t.perm))  ]
serum.t.perm <- serum.t.perm[, order(colMedians(serum.t.perm))  ]

serum.t.fdr <- serum.t.perm

for(i in 1 : nrow(serum.t.perm)){
  
  serum.t.fdr[i,] <- p.adjust(serum.t.perm[i,], method = "fdr")
}

serum.t.fdr <- serum.t.fdr[, order(colMedians(serum.t.fdr))  ]

write.csv(serum.t.perm, "Permutation_Result/serum_perm_p_vals.csv")
write.csv(serum.t.fdr, "Permutation_Result/serum_perm_p_vals_fdr.csv")
write.csv(serum.t.e.perm, "Permutation_Result/serum_perm_effect_sizes.csv")
write.csv(fold, "Permutation_Result/serum_perm_fold.csv")

################################################################################
## Make Figure 
stool.t.e.perm <- as.data.frame(stool.t.e.perm)

stool.t.fdr


serum.pick <- serum.t.e.perm[, colMedians(serum.t.fdr) < 0.05]
stool.t.e.perm[, colMedians(stool.t.fdr) < 0.05] # Choline is only significant metabolite after FDR

Choline <- stool.t.e.perm$Choline


mat <- cbind(serum.pick, Choline)

cols <- c("blue", "snow", "orangered")

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


layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
par(mar=c(12.1,4.1,4.1,2.1))

boxplot(mat,
        las = 2, 
        cex = .5,
        cex.axis=.8,
        pch = 21,
      #  ylim = c(-1.1, 1.1),
        ylab = "Effect size", 
        bg = "lightgrey",
        main = "METABOLITES 999 PERMUTATIONS",
        col = ff)

abline(h = 0, col = 2, lty = 2)
abline(v = 4.5, col = 1)
text(5, -1, "STOOL")
text(2.5, 1, "SERUM")

legend_image <- as.raster(matrix(colfunc(100), ncol=1))
plot(c(0,2),c(-1,1),type = 'n', axes = F,xlab = '', ylab = '',
     main = 'effect size', )

text(x=1.5, y = seq(-1,1,l=5), labels = seq(-1,1,l=5))
rasterImage(rev(legend_image), 0, -1, 1,1)
rect(0,-1,1,1)

segments(0,0,1,0, col = 2, lty = 2)


