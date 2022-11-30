# Sup_table_2
# Supplementary Table 2.
# Patients factors significantly correlated and/or associated with microbiome diversity.

source("R/human_functions.22.11.27.R")

map <- read.csv("data/Map_human_microbiome.22.11.21.csv")
row.names(map) <- map$sample.id

dummy <- grep("[-]", map$Collection.date)
keep.dummy <- as.Date(map$Collection.date[dummy], format = "%Y-%m-%d")
map$Collection.date <- as.Date(map$Collection.date, format = "%m/%d/%y")# dummy collection dates to exclude from analysis
map$Collection.date[dummy] <- keep.dummy

map.pat <- map[map$Disease.state == "Patient" , ]

ape.p <- NULL

age.dist2cent <- as.data.frame(matrix(nrow = 99, ncol = 4))
colnames(age.dist2cent) <- c("p value", "Spearman rho", "n", "S")
# S is the test statistic which is the sum of all squared rank differences

cor.perm <- function(var1 = map$eGFR,
                     var2 = map$richness,
                     donor = map$DonorID,
                     num.perm = 99,
                     norm = F){

  res <- as.data.frame(matrix(nrow =  num.perm, ncol = 4))
  colnames(res) <- c("p value", "Spearman rho", "n", "S")
  # S is the test statistic which is the sum of all squared rank differences

  pick <- rowSums(cbind(var1, var2))
  pick <- !is.na(pick)

  var1 <- var1[pick]
  var2 <- var2[pick]


  for(i in 1 : num.perm){

    # randomly select one sample from reseampled individuals
    pick2 <- sample(1:length(var1))

    var1 <- var1[pick2]
    var2 <- var2[pick2]
    donor <- donor[pick2]

    pick3 <- !duplicated(donor)

    if(norm ){
      res.i <- cor.test(var1[pick3], var2[pick3])
    }else{
      res.i <- cor.test(var1[pick3], var2[pick3], method = "spearman")
    }
    res$`p value`[i] <- res.i$p.value
    res$S[i] <- res.i$statistic
    res$`Spearman rho`[i] <- res.i$estimate
    res$n[i] <- length(var1[pick3])
  }

  print(res)

  return(sapply(res, mean))
}

# 1 AGE
cor.perm(map.pat$dist.cent, map.pat$age, map.pat$DonorID, 999)
## 2
cor.perm(map.pat$shared.proportion, map.pat$age, map.pat$DonorID, 999)
### 3
cor.perm(map.pat$dist2pair, map.pat$age, map.pat$DonorID, 999, norm = T)


#### 4 ORAL FOOD
inv.lme <- lmer(map.pat$invsimp ~ map.pat$Oral.foods + (1|map.pat$DonorID))
inv.null <- lmer(map.pat$invsimp ~ (1|map.pat$DonorID))
inv.res <- anova(inv.lme, inv.null)

inv.lme
inv.res

inv.sum <- c(inv.res$Df[2], inv.res$`Pr(>Chisq)`[2],  inv.res$Chisq[2],
             summary(inv.lme)$coefficients[2,1])
names(inv.sum) <- c("df", "p value", "Chi squared", "effect size")

##### 5
rich.lme <- lmer(map.pat$richness ~ map.pat$Oral.foods + (1|map.pat$DonorID))
rich.null <- lmer(map.pat$richness ~ (1|map.pat$DonorID))
rich.res <- anova(rich.lme, rich.null)

rich.lme
rich.res

rich.sum <- c(rich.res$Df[2], rich.res$`Pr(>Chisq)`[2],  rich.res$Chisq[2],
              summary(rich.lme)$coefficients[2,1])
names(rich.sum) <- c("df", "p value", "Chi squared", "effect size")

###### 6
pick <- !is.na(map.pat$shared.num)
shared.lme <- lmer(map.pat$shared.num[pick] ~ map.pat$Oral.foods[pick] + (1|map.pat$DonorID[pick]))
shared.null <- lmer(map.pat$shared.num[pick] ~ (1|map.pat$DonorID[pick]))
shared.res <- anova(shared.lme, shared.null)

shared.lme
shared.res

shared.sum <- c(shared.res$Df[2], shared.res$`Pr(>Chisq)`[2],  shared.res$Chisq[2],
              summary(shared.lme)$coefficients[2,1])
names(shared.sum) <- c("df", "p value", "Chi squared", "effect size")

####### 7 #########
fisher.test(table(map.pat$pam, map.pat$Oral.foods))

######## 8 ##########
fisher.test(table(map.pat$pam, map.pat$Antibiotic))

########### IMPORT METABOLOMICS data

meta.stool <- read.csv("data/HumanStoolUnpaired_biomassnorm_metalabels.csv", row.names = 1)
map.stool <- read.csv("data/Map_human_metabolome_stool.22.11.21.csv", row.names = 1)

meta.stool$Label <- c("dummy1", "dummy2")# needs multiple labels
pat.shared <- row.names(meta.stool)[row.names(meta.stool) %in% map.pat$sample.id]

pat.norm <- norm.human(metabolites = meta.stool, samples.pick = pat.shared)

m <- match(map.pat$sample.id, row.names(pat.norm))
# map.pat$sample.id == row.names(pat.norm)[m]
map.pat$stool.choline <- pat.norm$Choline[m]

### 9 #################
cor.perm(map.pat$invsimp, map.pat$stool.choline, map.pat$DonorID, 999, norm = T)


#### 10 #####################
cor.perm(map.pat$richness, map.pat$stool.choline, map.pat$DonorID, 999, norm = T)

##### 11 #####################
cor.perm(map.pat$shared.num, map.pat$stool.choline, map.pat$DonorID, 999, norm = T)


##### 12 ######################
plot(map.pat$shared.proportion, map.pat$Creatinine...35)
cor.perm(map.pat$shared.proportion, map.pat$Creatinine...35, map.pat$DonorID, 999, norm = F)

###### 13 ######################
plot(map.pat$dist2pair, map.pat$Creatinine...35)
cor.perm(map.pat$dist2pair, map.pat$Creatinine...35, map.pat$DonorID, 99, norm = F)


####### 14 ######################
cor.test(map.pat$eGFR, map.pat$richness)
cor.perm(log(map.pat$eGFR), map.pat$richness, map.pat$DonorID, 999, norm = T)

######## 15 ######################
cor.perm(log(map.pat$eGFR), map.pat$shared.num, map.pat$DonorID, 999, norm = T)

######### 16 ######################
cor.perm(log(map.pat$eGFR), map.pat$unique.num, map.pat$DonorID, 999, norm = T)

########## 17 ######################
cor.perm(map.pat$ALT, map.pat$shared.proportion, map.pat$DonorID, 999, norm = T)

########### 18 ######################
pick <- !is.na(map.pat$Funduplication) & !is.na(map.pat$stool.choline)
choline.lme <- lmer(map.pat$stool.choline[pick] ~ map.pat$Funduplication[pick] + (1|map.pat$DonorID[pick]))
choline.null <- lmer(map.pat$stool.choline[pick] ~ (1|map.pat$DonorID[pick]))
choline.res <- anova(choline.lme, choline.null)

choline.lme
choline.res

choline.sum <- c(choline.res$Df[2], choline.res$`Pr(>Chisq)`[2],  choline.res$Chisq[2],
                 summary(choline.lme)$coefficients[2,1])
names(choline.sum) <- c("df", "p value", "Chi squared", "effect size")

################### 19 ##################
meta.serum <- read.csv("data/HumanSerumUnpaired_biomassnorm_metalabels.csv", row.names = 1)
map.serum <- read.csv("data/Map_human_metabolome_serum.22.11.21.csv", row.names = 1)
map.serum$paired.stool
map.serum$BUN

map.serum.pat <- map.serum[map.serum$Disease.state == "Patient",]

meta.serum$Label <- c("dummy1", "dummy2")# needs multiple labels
serum.norm <- norm.human(metabolites = meta.serum, samples.pick = map.serum.pat$sample.id)
serum.norm <- norm.human(metabolites = meta.serum, samples.pick = map.serum$sample.id)


plot(meta.serum$Urea, log(map.serum$BUN))
cor.test(meta.serum$Urea, log(map.serum$BUN))

cor.perm(meta.serum$Urea, log(map.serum$BUN), map.pat$DonorID, 999, norm = F)



