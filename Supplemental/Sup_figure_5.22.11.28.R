getwd()
source("R/human_functions.22.11.27.R")

# makes figure 1 panel b
# Runs permutational T test
# this should take around 10 minutes

num.perm = 999

meta.stool <- read.csv("data/HumanStoolUnpaired_biomassnorm_metalabels.csv", row.names = 1)
meta.serum <- read.csv("data/HumanSerumUnpaired_biomassnorm_metalabels.csv", row.names = 1)

map.stool <- read.csv("data/Map_human_metabolome_stool.22.11.21.csv")
map.serum <- read.csv("data/Map_human_metabolome_serum.22.11.21.csv")

map.stool$Collection.date <- as.Date(map.stool$Collection.date, format = "%m/%d/%y")
map.serum$Collection.date <- as.Date(map.serum$Collection.date, format = "%m/%d/%y")

map.serum$col <- "#56B4E9" # assign colors be genotype
map.stool$col <- "#56B4E9" # assign colors be genotype

map.serum$col[map.serum$Disease.state == "Patient"] <- "#E69F00"
map.stool$col[map.stool$Disease.state == "Patient"] <- "#E69F00"

map.stool <- map.stool[match(row.names(meta.stool), map.stool$sample.id) , ]
map.serum <- map.serum[match(row.names(meta.serum), map.serum$sample.id) , ]

map.stool <- map.stool[!is.na(map.stool$Collection.date) , ]


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


#postscript("stool_p_t.test.eps", width = 16.0, height = 6.0)
par(mar=c(12.1,4.1,4.1,2.1))

boxplot(stool.t.fdr, las = 2, ylab = "p after FDR correction",
        main = paste("stool METABOLITES", num.perm, "PERMUTATIONS"))
abline(h = 0.05, col = 2)
#dev.off()
dir.create("Permutation_Result")

write.csv(stool.t.perm, "Permutation_Result/stool_perm_p_vals.csv")
write.csv(stool.t.fdr, "Permutation_Result/stool_perm_p_vals_fdr.csv")

write.csv(stool.t.e.perm, "Permutation_Result/stool_perm_effect_sizes.csv")
write.csv(stool.fold, "Permutation_Result/stool_perm_fold.csv")

boxplot(stool.t.e.perm)

aa <- colnames(stool.t.e.perm)
bb <- substr(aa, 1, 1)
cc <- substr(aa, 2, 2)
bb[bb == "X" & grepl("[0-9]", cc)] <- ""
dd <- paste0(bb, substr(aa, 2, 90))
ee <- gsub("[.]", "-", dd)

ee <- gsub("--", "-", ee)
ee <- gsub("-acid", " acid", ee)

colnames(stool.t.e.perm) <- ee

plot.effect(stool.t.e.perm, name = "stool_human_individual_metabolites",
            print = F)



########## compile source data file ###########################
permutation <- 1:999
source_data_5a <- cbind(permutation, stool.t.e.perm)

source_data_5a <- cbind("", "", source_data_5a)
source_data_5a[1,1] <- "figure 5b"
source_data_5a[1,2] <- "mean difference between paired samples"

source_data_5a <- rbind(source_data_5a, "", "")
write.csv(source_data_5a,"source_data_5a.csv", row.names = F)
#####################################################################

################################################################################
##                                                                            ##
##                                   SERUM                                    ##
##                                                                            ##
##                                                                            ##
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
meta.serum$X4.Pyridoxate <- NULL

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

  meta.serum.norm <- norm.human(samples.pick = samps.pick, metabolites = meta.serum)
  meta.serum.norm$Label <- NULL

  meta.i <- meta.serum.norm[match(map.i$sample.id, row.names(meta.serum.norm)) , ]

  for(i in 1 : ncol(meta.serum.norm)){

    s <-  meta.i[,i]

    res <- t.test(s[map.i$Disease.state == "Patient"],
                  s[map.i$Disease.state != "Patient"], paired = T)

    serum.t.perm[j,i] <- res$p.value
    serum.t.e.perm[j,i] <- res$estimate

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

fold <- cbind(comps, permutation, fold)

serum.t.e.perm <- serum.t.e.perm[, order(colMedians(serum.t.perm))  ]
serum.t.perm <- serum.t.perm[, order(colMedians(serum.t.perm))  ]

serum.t.fdr <- serum.t.perm

for(i in 1 : nrow(serum.t.perm)){

  serum.t.fdr[i,] <- p.adjust(serum.t.perm[i,], method = "fdr")
}

serum.t.fdr <- serum.t.fdr[, order(colMedians(serum.t.fdr))  ]



par(mar=c(12.1,4.1,4.1,2.1))

boxplot(serum.t.fdr, las = 2, ylab = "p after BH correction",
        main = paste("SERUM METABOLITES", num.perm, "PERMUTATIONS"))
abline(h = 0.05, col = 2)



write.csv(serum.t.perm, "Permutation_Result/serum_perm_p_vals.csv")
write.csv(serum.t.fdr, "Permutation_Result/serum_perm_p_vals_fdr.csv")


write.csv(serum.t.e.perm, "Permutation_Result/serum_perm_effect_sizes.csv")
write.csv(fold, "Permutation_Result/serum_perm_fold.csv")

aa <- colnames(serum.t.e.perm)
bb <- substr(aa, 1, 1)
cc <- substr(aa, 2, 2)
bb[bb == "X" & grepl("[0-9]", cc)] <- ""
dd <- paste0(bb, substr(aa, 2, 90))
ee <- gsub("[.]", "-", dd)

ee <- gsub("--", "-", ee)
ee <- gsub("-acid", " acid", ee)

colnames(serum.t.e.perm) <- ee

plot.effect(mat = serum.t.e.perm, name = "Serum_human_individual_metabolites",
            print = F)

########## compile source data file ###########################
permutation <- 1:999
source_data_5b <- cbind(permutation, serum.t.e.perm)

source_data_5b <- cbind("", "", source_data_5b)
source_data_5b[1,1] <- "figure 5a"
source_data_5b[1,2] <- "mean difference between paired samples"

write.csv(source_data_5b,"source_data_5b.csv", row.names = F)
#####################################################################
system("cat source_data_5a.csv source_data_5b.csv > source_data_5.csv")
system("rm source_data_5a.csv")
system("rm source_data_5b.csv")

