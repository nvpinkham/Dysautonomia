
getwd()
source("R/human_functions.22.06.18.R")

par(mfrow=c(1,5))

num.perm = 9

meta.stool <- read.csv("/Users/nickpinkham/Desktop/Nature2022_submission/data/human/HumanStoolUnpaired_biomassnorm_metalabels.csv", row.names = 1)
meta.serum <- read.csv("/Users/nickpinkham/Desktop/Nature2022_submission/data/human/HumanSerumUnpaired_biomassnorm_metalabels.csv", row.names = 1)

map.stool <- read.csv("/Users/nickpinkham/Desktop/Nature2022_submission/data/human/MasterStool12.15.21.csv")
map.serum <- read.csv("/Users/nickpinkham/Desktop/Nature2022_submission/data/human/MasterSerum12.15.21.csv")

map.stool$Collection.date <- as.Date(map.stool$Collection.date, format = "%m/%d/%y")
map.serum$Collection.date <- as.Date(map.serum$Collection.date, format = "%m/%d/%y")

map.serum$col <- "#56B4E9" # assign colors be genotype
map.stool$col <- "#56B4E9" # assign colors be genotype

map.serum$col[map.serum$Disease.state == "Patient"] <- "#E69F00"
map.stool$col[map.stool$Disease.state == "Patient"] <- "#E69F00"

map.stool <- map.stool[match(row.names(meta.stool), map.stool$sample.id) , ]
map.serum <- map.serum[match(row.names(meta.serum), map.serum$sample.id) , ]

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

stool.t.perm <-matrix(ncol = ncol(meta.stool) -1, nrow = num.perm)
colnames(stool.t.perm) <- colnames(meta.stool[-1])
stool.w.perm <- stool.diff <- stool.t.perm

map.stool$Collection.date <- as.Date(map.stool$Collection.date, "%m/%d/%y")

for(j in  1 : num.perm){

  print(j)

  map.i <- paired.map(map = map.stool, pats)

  ## DON'T RESAMPLE Individuals
  dubble.samps <- max(table(map.i$sample.id))
  while (   dubble.samps > 1 ) {
    print("resampleing because one or more samples were used
          more than once")

    map.i <- paired.map(map = map.stool, pats)
    dubble.samps <- max(table(map.i$sample.id))
  }

  meta.i<- norm.human(samples.pick = map.i$sample.id,
                      metabolites = meta.stool)


  meta.raw <- meta.stool[match( map.i$sample.id, row.names(meta.stool)) , ]

  row.names(meta.raw) == map.i$sample.id

  meta.i$Label <- NULL
  meta.raw$Label <- NULL

  for(i in 1 : ncol(meta.i)){

    s <-  meta.i[,i]
    s.raw <- meta.raw[,i]

    colnames(meta.i)[i]

    p <- log(s.raw[map.i$Disease.state == "Patient"])
    r <- log(s.raw[map.i$Disease.state == "Relative"])

    res <- t.test(s[map.i$Disease.state == "Patient"],
                  s[map.i$Disease.state != "Patient"], paired = T)

    stool.t.perm[j,i] <- res$p.value

    #stool.diff[j,i] <- mean(s[map.i$Disease.state == "Patient"] - s[map.i$Disease.state != "Patient"] / s[map.i$Disease.state != "Patient"]) * 100
    stool.diff[j,i] <- mean( (p - r) / r ) * 100

    #stool.diff[j,i] <- mean(log(s.raw[map.i$Disease.state == "Patient"])/
    #                         log(s.raw[map.i$Disease.state != "Patient"])) * 100

  }
}

colnames(stool.t.perm)[colMeans(stool.t.perm) < .05]

all(map.stool$sample.id == row.names(meta.stool))

stool.t.fdr <- stool.t.perm

for(i in 1 : nrow(stool.t.perm)){

  stool.t.fdr[i,] <- p.adjust(stool.t.perm[i,], method = "fdr")

}

sort(colMeans(stool.t.fdr))[1:5]

qqnorm(meta.stool$Choline)
qqnorm(log(meta.stool$Choline))

shapiro.test(meta.stool$Choline)
shapiro.test(log(meta.stool$Choline))

stool.diff <- stool.diff[ , order(colMeans(stool.t.fdr))]
stool.t.fdr <- stool.t.fdr[ , order(colMeans(stool.t.fdr))]

par(mfrow=c(1,5))


# 58 pairs
a <- log(meta.stool$Choline[map.stool$Disease.state == "Patient"])
b <- log(meta.stool$Choline[map.stool$Disease.state != "Patient"])

mean(a)
mean(b)

(mean(b) - mean(a)) / mean(a) * 100

table(map.serum$Disease.state)

p <- meta.i$Choline[map.i$Disease.state == "Patient"]
r <- meta.i$Choline[map.i$Disease.state == "Relative"]

p <- meta.raw$Choline[map.i$Disease.state == "Patient"]
r <- meta.raw$Choline[map.i$Disease.state == "Relative"]


stool.nrom <- norm.human(map.stool$sample.id, meta.stool)

par(mar=c(4.1, 4.1, 5.1, 2.1))

paired.violin.beta(map = map.stool, var = log(meta.stool$Choline))
#paired.violin.beta(map = map.stool, var = stool.nrom$Choline)
title(paste("Stool Choline", "\nmean t test after FDR\n p value =",
            round(mean(stool.t.fdr[,1]), 5),
            "\nmean % change =", abs(round(mean(stool.diff[,1]), 2))))
title(ylab = "natural log")

par(mar=c(4.1, 3.1, 5.1, 2.1))

########### SERUM



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
norm.pat <- norm.rel <- serum.w.perm <- serum.diff <- serum.t.perm

map.serum$Collection.date <- as.Date(map.serum$Collection.date, "%m/%d/%y")

for(j in  1 : num.perm){

  print(j)

  map.i <- paired.map(map = map.serum, pats)

  ## DON'T RESAMPLE Individuals
  dubble.samps <- max(table(map.i$sample.id))
  while (   dubble.samps > 1 ) {
    print("resampleing because one or more samples were used
          more than once")

    map.i <- paired.map(map = map.serum, pats)
    dubble.samps <- max(table(map.i$sample.id))
  }

  meta.i<- norm.human(samples.pick = map.i$sample.id,
                      metabolites = meta.serum)


  meta.raw <- meta.serum[match( map.i$sample.id, row.names(meta.serum)) , ]

  row.names(meta.raw) == map.i$sample.id

  meta.i$Label <- NULL
  meta.raw$Label <- NULL

  for(i in 1 : ncol(meta.i)){

    s <-  meta.i[,i]
    s.raw <- meta.raw[,i]

    colnames(meta.i)[i]

    p <- log(s.raw[map.i$Disease.state == "Patient"])
    r <- log(s.raw[map.i$Disease.state == "Relative"])

    norm.pat[j,i] <- shapiro.test(s[map.i$Disease.state == "Patient"])$p.value
    norm.rel[j,i] <- shapiro.test(s[map.i$Disease.state != "Patient"])$p.value

    res <- t.test(s[map.i$Disease.state == "Patient"],
                  s[map.i$Disease.state != "Patient"], paired = T)

    serum.t.perm[j,i] <- res$p.value

    serum.diff[j,i] <- mean( (p - r) / r ) * 100

  }
}

colnames(serum.t.perm)[colMeans(serum.t.perm) < .05]

all(map.serum$sample.id == row.names(meta.serum))

colnames(serum.t.perm)[colMeans(norm.pat) < .05]
colnames(serum.t.perm)[colMeans(norm.rel) < .05]

serum.t.fdr <- serum.t.perm
serum.w.fdr <- serum.w.perm

for(i in 1 : nrow(serum.t.perm)){

  serum.t.fdr[i,] <- p.adjust(serum.t.perm[i,], method = "fdr")
  serum.w.fdr[i,] <- p.adjust(serum.w.perm[i,], method = "fdr")

}

sort(colMeans(serum.t.fdr))[1:5]

serum.diff <- serum.diff[ , order(colMeans(serum.t.fdr))]
serum.t.fdr <- serum.t.fdr[ , order(colMeans(serum.t.fdr))]

table(map.serum$Disease.state)

paired.violin.beta(map = map.serum, var = log(meta.serum$X4.Pyridoxate))

title(paste("Serum X4.Pyridoxate",
            "\nmean t test after FDR\n p value =",
            round(mean(serum.t.fdr[,1]), 5),
            "\nmean% change =", abs(round(mean(serum.diff[,1]), 2))))


paired.violin.beta(map = map.serum, var = log(meta.serum$Xanthine))

title(paste("Serum Xanthine", "\nmean t test after FDR\n p value =", round(mean(serum.t.fdr[,2]), 5),
            "\nmean% change =", abs(round(mean(serum.diff[,2]), 2))))


paired.violin.beta(map = map.serum, var = log(meta.serum$Urea))

title(paste("Serum Urea", "\nmean t test after FDR\n p value =", round(mean(serum.t.fdr[,3]), 5),
            "\nmean% change =", abs(round(mean(serum.diff[,3]), 2))))


paired.violin.beta(map = map.serum, var = log(meta.serum$Methanol))

title(paste("Serum Methanol", "\nmean t test after FDR\n p value =", round(mean(serum.t.fdr[,4]), 5),
            "\nmean% change =", abs(round(mean(serum.diff[,4]), 2))))

