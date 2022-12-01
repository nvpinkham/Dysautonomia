
source("R/human_functions.22.11.27.R")

# makes sup figure 5
# Runs permutational T test
# this should take around 20 minutes

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
stool.t.hi.perm <- stool.t.lo.perm <- stool.t.t.perm <- stool.t.e.perm <- stool.t.perm.no <- stool.t.perm

fold <-  matrix(ncol = ncol(meta.stool) - 1, nrow = num.perm * length(pats))

d <- length(pats)
permutation <- NULL
comps <- NULL
comps.list <- list()

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

    stool.t.e.perm[j,i] <- res$estimate
    stool.t.t.perm[j,i] <- res$statistic
    stool.t.lo.perm[j,i] <- res$conf.int[1]
    stool.t.hi.perm[j,i] <- res$conf.int[2]

    here <-  (((j - 1) * d) + 1) : (j * d)

    fold[ here , i] <- s[map.i$Disease.state == "Patient"] /s[map.i$Disease.state != "Patient"]

    permutation[here] <- rep(j, d)

    comp <- paste(map.i$sample.id[map.i$Disease.state == "Patient"],
                  "to",
                  map.i$sample.id[map.i$Disease.state != "Patient"])


    comps[here] <- comp

    comps.list[[j]] <- comp
  }
}

stool.fold <- cbind(comps, permutation, fold)

stool.t.e.perm <- stool.t.e.perm[, order(colMedians(stool.t.perm))  ]
stool.t.perm <- stool.t.perm[, order(colMedians(stool.t.perm))  ]

stool.t.fdr <- stool.t.perm

for(i in 1 : nrow(stool.t.perm)){

  stool.t.fdr[i,] <- p.adjust(stool.t.perm[i,], method = "fdr")
}



`mean p value` <- apply(stool.t.perm, 2, mean)
`mean p after FDR` <- apply(stool.t.fdr, 2, mean)
`mean effect size` <- apply(stool.t.e.perm, 2, mean)
`mean CI 95% lo` <- apply(stool.t.lo.perm, 2, mean)
`mean CI 95% high` <- apply(stool.t.hi.perm, 2, mean)
`mean df` <- res$parameter
`mean t stat` <- apply(stool.t.t.perm, 2, mean)

stool.res <- cbind(`mean p value`,
                   `mean p after FDR`,
                   `mean effect size` ,
                   `mean CI 95% lo` ,
                   `mean CI 95% high` ,
                   `mean df`,
                   `mean t stat`)

########## compile source data file ###########################
permutation <- 1 : num.perm
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
serum.t.hi.perm <- serum.t.lo.perm <- serum.t.t.perm <- serum.t.e.perm <- serum.t.perm.no <- serum.t.perm

#eff <- array(dim = c(ncol(meta.serum), length(pats), num.perm))
fold <-  matrix(ncol = ncol(meta.serum) - 1, nrow = num.perm * length(pats))
dim(fold)

d <- length(pats)
permutation <- NULL
comps <- NULL
comps.list <- list()

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

    serum.t.e.perm[j,i] <- res$estimate
    serum.t.t.perm[j,i] <- res$statistic
    serum.t.lo.perm[j,i] <- res$conf.int[1]
    serum.t.hi.perm[j,i] <- res$conf.int[2]

    here <-  (((j - 1) * d) + 1) : (j * d)

    fold[ here , i] <- s[map.i$Disease.state == "Patient"] /s[map.i$Disease.state != "Patient"]

    permutation[here] <- rep(j, d)

    comp <- paste(map.i$sample.id[map.i$Disease.state == "Patient"],
                  "to",
                  map.i$sample.id[map.i$Disease.state != "Patient"])


    comps[here] <- comp

    comps.list[[j]] <- comp
  }
}

fold <- cbind(comps, permutation, fold)

serum.t.e.perm <- serum.t.e.perm[, order(colMedians(serum.t.perm))  ]
serum.t.perm <- serum.t.perm[, order(colMedians(serum.t.perm))  ]

serum.t.fdr <- serum.t.perm

for(i in 1 : nrow(serum.t.perm)){

  serum.t.fdr[i,] <- p.adjust(serum.t.perm[i,], method = "fdr")
}

serum.t.fdr <- serum.t.fdr[, order(colMedians(serum.t.fdr))  ]

`mean p value` <- apply(serum.t.perm, 2, mean)
`mean p after FDR` <- apply(serum.t.fdr, 2, mean)
`mean effect size` <- apply(serum.t.e.perm, 2, mean)
`mean CI 95% lo` <- apply(serum.t.lo.perm, 2, mean)
`mean CI 95% high` <- apply(serum.t.hi.perm, 2, mean)
`mean df` <- res$parameter
`mean t stat` <- apply(serum.t.t.perm, 2, mean)

serum.res <- cbind(`mean p value`,
                   `mean p after FDR`,
                   `mean effect size` ,
                   `mean CI 95% lo` ,
                   `mean CI 95% high` ,
                   `mean df`,
                   `mean t stat`)

choline <- stool.res[stool.res[,2] < 0.05 , ] # really good trick
sup_table_3 <- rbind(choline,
                     serum.res[serum.res[,2] < 0.05 , ])

write.csv(sup_table_3, "Supplemental_table_3.csv")



