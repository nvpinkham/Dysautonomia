
getwd()
source("R/human_functions.22.11.21.R")

par(mfrow=c(1,4))

num.perm = 999

meta.stool <- read.csv("data/HumanStoolUnpaired_biomassnorm_metalabels.csv", row.names = 1)
meta.serum <- read.csv("data/HumanSerumUnpaired_biomassnorm_metalabels.csv", row.names = 1)

map.stool <- read.csv("data/Map_human_metabolome_stool.csv", row.names = 1)
map.serum <- read.csv("data/Map_human_metabolome_serum.csv", row.names = 1)

map.stool$Collection.date <- as.Date(map.stool$Collection.date, format = "%Y-%m-%d")# excel messes up dates by careful
map.serum$Collection.date <- as.Date(map.serum$Collection.date, format = "%Y-%m-%d")

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

# savoutput of t.test for each permuation
stool.t.perm <-matrix(ncol = ncol(meta.stool) -1, nrow = num.perm)
colnames(stool.t.perm) <- colnames(meta.stool[-1])

stool.t.perm.df <- stool.t.perm
stool.t.perm.effect_size <- stool.t.perm
stool.t.perm.t.stat <- stool.t.perm
stool.t.perm.conf_int_lo <- stool.t.perm
stool.t.perm.conf_int_hi <- stool.t.perm


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

  meta.i <- meta.i[match( map.i$sample.id, row.names(meta.i)) , ]

  all(row.names(meta.i) == map.i$sample.id)

  meta.i$Label <- NULL

  for(i in 1 : ncol(meta.i)){

    s <-  meta.i[,i]

    res <- t.test(s[map.i$Disease.state == "Patient"],
                  s[map.i$Disease.state != "Patient"], paired = T)

    stool.t.perm[j,i] <- res$p.value
    stool.t.perm.df[j,i] <- res$parameter
    stool.t.perm.effect_size[j,i] <- res$estimate
    stool.t.perm.t.stat[j,i] <- res$statistic
    stool.t.perm.conf_int_lo[j,i] <- res$conf.int[1]
    stool.t.perm.conf_int_hi[j,i] <- res$conf.int[2]
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


stool.t.fdr <- stool.t.fdr[ , order(colMeans(stool.t.fdr))]

par(mfrow=c(1,4))


par(mar=c(4.1, 4.1, 5.1, 2.1))

paired.violin.beta(map = map.stool, var = log(meta.stool$Choline))
#paired.violin.beta(map = map.stool, var = stool.nrom$Choline)
title(paste("Stool Choline", "\nmean t test after FDR\n p value =",
            round(mean(stool.t.fdr[,1]), 5)))
title(ylab = "natural log")

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

serum.t.perm.df <- serum.t.perm
serum.t.perm.effect_size <- serum.t.perm
serum.t.perm.t.stat <- serum.t.perm
serum.t.perm.conf_int_lo <- serum.t.perm
serum.t.perm.conf_int_hi <- serum.t.perm


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

  meta.i <- meta.i[match( map.i$sample.id, row.names(meta.i)) , ]

  all(row.names(meta.i) == map.i$sample.id)

  meta.i$Label <- NULL

  for(i in 1 : ncol(meta.i)){

    s <-  meta.i[,i]

    norm.pat[j,i] <- shapiro.test(s[map.i$Disease.state == "Patient"])$p.value
    norm.rel[j,i] <- shapiro.test(s[map.i$Disease.state != "Patient"])$p.value

    res <- t.test(s[map.i$Disease.state == "Patient"],
                  s[map.i$Disease.state != "Patient"], paired = T)

    serum.t.perm[j,i] <- res$p.value
    serum.t.perm.df[j,i] <- res$parameter
    serum.t.perm.effect_size[j,i] <- res$estimate
    serum.t.perm.t.stat[j,i] <- res$statistic
    serum.t.perm.conf_int_lo[j,i] <- res$conf.int[1]
    serum.t.perm.conf_int_hi[j,i] <- res$conf.int[2]

    # serum.diff[j,i] <- mean( (p - r) / r ) * 100
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

paired.violin.beta(map = map.serum, var = log(meta.serum$Xanthine))

title(paste("Serum Xanthine", "\nmean t test after FDR\n p value =", round(mean(serum.t.fdr[,2]), 5),
            "\nmean% change =", abs(round(mean(serum.diff[,2]), 2))))


paired.violin.beta(map = map.serum, var = log(meta.serum$Urea))

title(paste("Serum Urea", "\nmean t test after FDR\n p value =", round(mean(serum.t.fdr[,3]), 5),
            "\nmean% change =", abs(round(mean(serum.diff[,3]), 2))))


paired.violin.beta(map = map.serum, var = log(meta.serum$Methanol))

title(paste("Serum Methanol", "\nmean t test after FDR\n p value =", round(mean(serum.t.fdr[,4]), 5),
            "\nmean% change =", abs(round(mean(serum.diff[,4]), 2))))

########################################

source_data_2b <- cbind(map.stool$sample.id,
                        map.stool$DonorID,
                        map.stool$Family.ID,
                        map.stool$Disease.state)

colnames(source_data_2b) <- c("sample.id", "DonorID", "Family.ID", "Disease.state")

source_data_2b <- cbind("", "", source_data_2b, "","",
                        log(meta.stool$Choline))

source_data_2b[1,1] <- "Figure 2a"
source_data_2b[1,2] <- "Stool"
source_data_2b[1,8] <- "natural log stool metablite concentration"

colnames(source_data_2b)[8] <- "Choline"

source_data_2b <- rbind(source_data_2b, "", "")
write.csv(source_data_2b,"source_data_2b1.csv", row.names = F)


source_data_2b <- cbind(map.serum$sample.id,
                        map.serum$DonorID,
                        map.serum$Family.ID,
                        map.serum$Disease.state)

colnames(source_data_2b) <- c("sample.id", "DonorID", "Family.ID", "Disease.state")

source_data_2b <- cbind("", "", source_data_2b, "","",
                        log(meta.serum$Xanthine),
                        log(meta.serum$Urea),
                        log(meta.serum$Methanol))

source_data_2b[1,2] <- "Serum"
source_data_2b[1,8] <- "natural log serum metablite concentration"

colnames(source_data_2b)[8:10] <- c("Xanthine", "Urea", "Methanol")

write.csv(source_data_2b,"source_data_2b2.csv", row.names = F)

system("cat source_data_2b1.csv source_data_2b2.csv > source_data_2b.csv")





