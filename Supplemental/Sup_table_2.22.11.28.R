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

# 1 AGE

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

cor.perm(map.pat$dist.cent, map.pat$age, map.pat$DonorID, 999)
cor.perm(map.pat$shared.proportion, map.pat$age, map.pat$DonorID, 999)
cor.perm(map.pat$dist2pair, map.pat$age, map.pat$DonorID, 999, norm = T)


#2 ORAL FOOD
inv.lme <- lmer(map.pat$invsimp ~ map.pat$Oral.foods + (1|map.pat$DonorID))
inv.null <- lmer(map.pat$invsimp ~ (1|map.pat$DonorID))
inv.lme

anova(inv.lme, inv.null)

rich.lme <- lmer(map.pat$richness ~ map.pat$Oral.foods + (1|map.pat$DonorID))
rich.null <- lmer(map.pat$richness ~ (1|map.pat$DonorID))

anova(rich.lme, rich.null)
rich.lme

pick <- !is.na(map.pat$shared.num)
shared.lme <- lmer(map.pat$shared.num[pick] ~ map.pat$Oral.foods[pick] + (1|map.pat$DonorID[pick]))
shared.null <- lmer(map.pat$shared.num[pick] ~ (1|map.pat$DonorID[pick]))

anova(shared.lme, shared.null)
shared.lme
aa <- summary(shared.lme)
diff(aa$coefficients[,1])


fisher.test(map.pat$pam, map.pat$Oral.foods)

fisher.test(map.pat$pam, map.pat$Antibiotic)
####

cor.perm(map.pat$eGFR, map.pat$richness, map.pat$DonorID, 99, norm = T)
cor.perm(map$eGFR, map$richness, map$DonorID, 99, norm = T)

table(!is.na(map$eGFR))
table(!is.na(map$eGFR[map$Disease.state == "Patient"]))

#####
table(is.na( ! map$eGFR))

table(map.pat$DonorID)

map.pat1 <- map.pat[sample(1:nrow(map.pat)), ]
map.unique <- map.pat1[!duplicated(map.pat1$DonorID) , ]
sort(map.pat$Collection.date)
sort(map.unique$Collection.date)

cor.test(map.unique$richness, map.unique$eGFR)
cor.test(map.pat$shared.num, map.pat$eGFR)
cor.test(map$unique.num, map$eGFR)

pick <- !is.na(map.pat$shared.num) & !is.na(map.pat$eGFR)
shared.lme <- lmer(map.pat$shared.num[pick] ~ map.pat$eGFR[pick] + (1|map.pat$DonorID[pick]))
shared.null <- lmer(map.pat$shared.num[pick] ~ (1|map.pat$DonorID[pick]))
anova(shared.lme, shared.null)


cor.test(map$richness, map$eGFR)
cor.test(map$richness, map$eGFR)



fish.res <- fisher.test(table(map$enterotype[map$Disease.state == "Patient"],
                              map$G.tube[map$Disease.state == "Patient"]))

