source("R/human_functions.22.11.27.R")

# makes figure 1 panel c
# Runs permutational T test
# this should take around 10 minutes

num.perm = 999

map <- read.csv("data/Map_human_microbiome.22.11.21.csv")
row.names(map) <- map$sample.id

dummy <- grep("[-]", map$Collection.date)
keep.dummy <- as.Date(map$Collection.date[dummy], format = "%Y-%m-%d")
map$Collection.date <- as.Date(map$Collection.date, format = "%m/%d/%y")# dummy collection dates to exclude from analysis
map$Collection.date[dummy] <- keep.dummy

otu <- read.csv("data/FD_OTU_human21.csv", row.names = 1)

reconcile(map, otu)

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

otu <- otu[filter,]
map <- map[filter,]

reconcile(map, otu)

table(map$Disease.state)
################################################################################

pats <-  unique(map$DonorID[map$Disease.state == "Patient"])

alpha.t <- as.data.frame(matrix(nrow = num.perm, ncol = 6))

rich.t <- beta.t <- alpha.t

colnames(alpha.t) <- c("P_val", "df", "effect_size", "t.stat", "conf_int_lo",  "conf_int_hi")
colnames(beta.t) <- c("P_val", "df", "effect_size", "t.stat", "conf_int_lo",  "conf_int_hi")
colnames(rich.t) <- c("P_val", "df", "effect_size", "t.stat", "conf_int_lo",  "conf_int_hi")

families <- unique(map.i$Family.ID)

if(!exists("otu.ps")){
  # only run if result is not in the environment

  print("new permutation analysis")
  otu.effect <- otus <- otu.ps <- list()

  otus.df <- as.data.frame(matrix(ncol = ncol(otu),
                                  nrow = num.perm))

  colnames(otus.df) <- colnames(otu)
  otus.ts <- otus.hi <- otus.lo <- otus.ps.bh <- otus.ps <- otus.ef <- otus.df

  otus <- otu.ps <- otu.effect <- list()
  proportion.p <- proportion.ef <- unique.p <- unique.ef <- NULL
  p.mean <- r.mean <- NULL
  pat.disp <- rel.disp <- NULL

  donors <- unique(map$DonorID)
  donors.uni <- matrix(ncol = length(donors), nrow = num.perm)
  colnames(donors.uni) <- sort(donors)

  for(i in  1 : num.perm){

    print(i)
    map.i <- paired.map(map, pats)

    res <- t.test(map.i$invsimp ~ map.i$Disease.state, paired = T)

    alpha.t$P_val[i] <- res$p.value
    alpha.t$df[i] <- res$parameter
    alpha.t$t.stat[i] <- res$statistic
    alpha.t$effect_size[i] <- res$estimate
    alpha.t$conf_int_lo[i] <- res$conf.int[1]
    alpha.t$conf_int_hi[i] <- res$conf.int[2]

    res <- t.test(map.i$richness ~ map.i$Disease.state, paired = T)

    rich.t$P_val[i] <- res$p.value
    rich.t$df[i] <- res$parameter
    rich.t$t.stat[i] <- res$statistic
    rich.t$effect_size[i] <- res$estimate
    rich.t$conf_int_lo[i] <- res$conf.int[1]
    rich.t$conf_int_hi[i] <- res$conf.int[2]

    ## explore beta diversity
    otu.i <- otu[match(map.i$sample.id, row.names(otu)) , ]

    otu.pick <- otu.i[ , colSums(otu.i > 0) > 0]

    # b <- matrix(nrow = length( families), ncol =  ncol (otu.pick))

    # leverage the paired nature of the experiment
    # paired t test of each OTU

    ################################
    res <- sapply(otu.pick,  function(otu.pick) t.test(otu.pick ~ map.i$Disease.state, paired = T))

    dim(res)

    ps <- res [3,]
    ps <- sapply(  ps, as.numeric)

    ps.bh <- p.adjust(ps, method = "BH")

    row.names(res)

    otus.ts[i, match(names(ps), colnames(otus.df)) ] <- res[1,]
    otus.df[i, match(names(ps), colnames(otus.df)) ] <- res[2,]
    otus.ps[i, match(names(ps), colnames(otus.df)) ] <- res[3 ,]
    otus.lo[i, match(names(ps), colnames(otus.df)) ] <- sapply(res[4 ,], "[[", 1)
    otus.hi[i, match(names(ps), colnames(otus.df)) ] <- sapply(res[4 ,], "[[", 2)
    otus.ef[i, match(names(ps), colnames(otus.df)) ] <- res[5 ,]

    ##########
    b.res <- betadisper(vegdist(otu.i), map.i$Disease.state)

    pat.disp[[i]] <- b.res$distances[ map.i$Disease.state == "Patient"]
    rel.disp[[i]] <- b.res$distances[ map.i$Disease.state != "Patient"]

    res <- t.test(b.res$distances ~ b.res$group, paired = T)# distance to group centroid

    beta.t$P_val[i] <- res$p.value
    beta.t$df[i] <- res$parameter
    beta.t$effect_size[i] <- res$estimate
    beta.t$t.stat[i] <- res$statistic
    beta.t$conf_int_lo[i] <- res$conf.int[1]
    beta.t$conf_int_hi[i] <- res$conf.int[2]

  }
}

########### OTUS mentioned in the text #############
colnames(otu)

otu.summaries <- cbind(colMeans(otus.ps, na.rm = T),
                       colMeans(otus.ef, na.rm = T),
                       colMeans(otus.ts, na.rm = T),
                       colMeans(otus.df, na.rm = T),
                       colMeans(otus.lo, na.rm = T),
                       colMeans(otus.hi, na.rm = T))

colnames(otu.summaries) <- c("mean ps value",
                             "mean effect size",
                             "mean t stat",
                             "df",
                             "lo 95% CI",
                             "hi 95% CI")

otu.sum.pick <- otu.summaries[c(4, 9),]

################################################################################

dir.create("Microbiome_Permution_Analysis")
saveRDS(alpha.t, "Microbiome_Permution_Analysis/alpha.t.rds")
saveRDS(rich.t, "Microbiome_Permution_Analysis/rich.t.rds")
saveRDS(beta.t, "Microbiome_Permution_Analysis/beta.t.rds")
write.csv(otu.summaries, "Microbiome_Permution_Analysis/all_OTU_t_tests.csv")


################################################################################
#####            RICHNESS                  ###############

rich.t <- readRDS("Microbiome_Permution_Analysis/rich.t.rds")

paired.violin.beta(map, var = map$richness)
title(ylab ="Richness")

title(paste("Paired T test mean p value =",  round(mean(rich.t$P_val, na.rm = T), 5),
            "\nmean effect size =", round(mean(rich.t$effect_size, na.rm = T), 5),
            "\nmean 95% CI =", round(mean(rich.t$conf_int_lo, na.rm = T), 2), "to",
            round(mean(rich.t$conf_int_hi, na.rm = T), 2),
            "\ndf =", round(mean(rich.t$df, na.rm = T), 5)))


##################################################################

source_data.1c <- cbind(map$sample.id,  map$DonorID, map$Family.ID, map$Disease.state, map$richness)
colnames(source_data.1c) <- c("sample.id",  "DonorID", "Family.ID", "Disease.state", "Richness")

source_data.1c <- cbind("", source_data.1c)
source_data.1c[1,1] <- "Figure 1c"

write.csv(source_data.1c, "source_data/source_data_1c.csv")
rich.t.means <- apply(rich.t, 2, mean)

names(rich.t.means) <- paste("mean over", num.perm, "permutations", names(rich.t.means))

rich.t.means <- capture.output(rich.t.means)

writeLines(rich.t.means, "Statistical_summaries/Fig_1c_tests.txt")

