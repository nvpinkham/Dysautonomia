
source("R/human_functions.22.11.27.R")

par(mfrow=c(1,4))

num.perm = 999

meta.stool <- read.csv("data/HumanStoolUnpaired_biomassnorm_metalabels.csv", row.names = 1)
meta.serum <- read.csv("data/HumanSerumUnpaired_biomassnorm_metalabels.csv", row.names = 1)

map.stool <- read.csv("data/Map_human_metabolome_stool.22.11.21.csv", row.names = 1)
map.serum <- read.csv("data/Map_human_metabolome_serum.22.11.21.csv", row.names = 1)

map.stool$Collection.date <- as.Date(map.stool$Collection.date, format = "%m/%d/%y")# excel messes up dates by careful
map.serum$Collection.date <- as.Date(map.serum$Collection.date, format = "%m/%d/%y")

map.stool$sample.id <- row.names(map.stool)

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

ranking <- order(colMeans(stool.t.fdr))
stool.t.fdr <- stool.t.fdr[ , ranking]
stool.t.perm.effect_size <- stool.t.perm.effect_size[ , ranking ]
stool.t.perm.conf_int_hi <- stool.t.perm.conf_int_hi[ , ranking ]
stool.t.perm.conf_int_lo <- stool.t.perm.conf_int_lo[ , ranking ]


par(mfrow=c(1,4))


par(mar=c(4.1, 4.1, 9.1, 2.1))

paired.violin.beta(map = map.stool, var = log(meta.stool$Choline))
#paired.violin.beta(map = map.stool, var = stool.nrom$Choline)
title(paste("Stool Choline",
            "\n t test after FDR\n mean p value =",
            round(mean(stool.t.fdr[,1]), 5),
            "\nmean effect size =", round(mean(stool.t.perm.effect_size[,1]), 5),
            "\nmean 95% CI =", round(mean(stool.t.perm.conf_int_lo[,1]), 3),
            round(mean(stool.t.perm.conf_int_hi[,1]), 3)))

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


ranking <- order(colMeans(serum.t.fdr))
serum.t.fdr <- serum.t.fdr[ , ranking]
serum.t.perm.effect_size <- serum.t.perm.effect_size[ , ranking ]
serum.t.perm.conf_int_hi <- serum.t.perm.conf_int_hi[ , ranking ]
serum.t.perm.conf_int_lo <- serum.t.perm.conf_int_lo[ , ranking ]


paired.violin.beta(map = map.serum, var = log(meta.serum$Xanthine))

title(paste("Serum Xanthine",
            "\n t test after FDR\n mean p value =",
            round(mean(serum.t.fdr[,1]), 5),
            "\nmean effect size =", round(mean(serum.t.perm.effect_size[,1]), 5),
            "\nmean 95% CI =", round(mean(serum.t.perm.conf_int_lo[,1]), 3),
            round(mean(serum.t.perm.conf_int_hi[,1]), 3)))

paired.violin.beta(map = map.serum, var = log(meta.serum$Urea))

title(paste("Serum Urea",
            "\n t test after FDR\n mean p value =",
            round(mean(serum.t.fdr[,2]), 5),
            "\nmean effect size =", round(mean(serum.t.perm.effect_size[,2]), 5),
            "\nmean 95% CI =", round(mean(serum.t.perm.conf_int_lo[,2]), 3),
            round(mean(serum.t.perm.conf_int_hi[,2]), 3)))

paired.violin.beta(map = map.serum, var = log(meta.serum$Methanol))

title(paste("Serum Methanol",
            "\n t test after FDR\n mean p value =",
            round(mean(serum.t.fdr[,3]), 5),
            "\nmean effect size =", round(mean(serum.t.perm.effect_size[,3]), 5),
            "\nmean 95% CI =", round(mean(serum.t.perm.conf_int_lo[,3]), 3),
            round(mean(serum.t.perm.conf_int_hi[,3]), 3)))

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
write.csv(source_data_2b,"source_data/source_data_2b1.csv", row.names = F)


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

write.csv(source_data_2b,"source_data/source_data_2b2.csv", row.names = F)

system("cat source_data/source_data_2b1.csv source_data/source_data_2b2.csv > source_data/source_data_2b.csv")
system("rm source_data/source_data_2b1.csv")
system("rm source_data/source_data_2b2.csv")


#######

stool.res <- cbind(apply(stool.t.fdr, 2, mean),
                   apply(stool.t.perm.effect_size, 2, mean),
                   apply(stool.t.perm.df, 2, mean),
                   paste(apply(stool.t.perm.conf_int_lo, 2, mean),
                         "to",
                         apply(stool.t.perm.conf_int_hi, 2, mean)))

colnames(stool.res) <- c("Stool mean p value after FDR",
                         "mean effect size",
                         "df",
                         "mean 95% CI")

write.table(stool.res, "Statistical_summaries/Fig_2b_stool_tests.txt", quote = F)



serum.res <- cbind(apply(serum.t.fdr, 2, mean),
                   apply(serum.t.perm.effect_size, 2, mean),
                   apply(serum.t.perm.df, 2, mean),
                   paste(apply(serum.t.perm.conf_int_lo, 2, mean),
                         "to",
                         apply(serum.t.perm.conf_int_hi, 2, mean)))

colnames(serum.res) <- c("Serum mean p value after FDR",
                         "mean effect size",
                         "df",
                         "mean 95% CI")

write.table(serum.res, "Statistical_summaries/Fig_2b_serum_tests.txt", quote = F)


system("cat Statistical_summaries/Fig_2b_stool_tests.txt Statistical_summaries/Fig_2b_serum_tests.txt > Statistical_summaries/Fig_2b_tests.txt")
system("rm Statistical_summaries/Fig_2b_stool_tests.txt")
system("rm Statistical_summaries/Fig_2b_serum_tests.txt")



