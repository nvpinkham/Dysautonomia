
source("R/mice_functions.22.12.01.R")

map <- read.csv("data/Map_CoHoused_mice.8.16.csv")
map$discription <- paste0(map$Genotype, "_", map$Treatment.type)
map$ID.Gen.Tmt.Age <- paste0("X", map$ID.Gen.Tmt.Age)
map$DOB <- as.Date(map$DOB , "%m/%d/%y")

map$X16SFile.ID[map$X16SFile.ID == "" ] <- NA

map$col[map$Treatment.type =="Separate" & map$Genotype == "Mutant" ] <- '#E69F00'
map$col[map$Treatment.type =="Separate" & map$Genotype == "Control"] <- '#56B4E9'
map$col[map$Treatment.type =="Cohoused" & map$Genotype == "Mutant"] <-  '#009E73'
map$col[map$Treatment.type =="Cohoused" & map$Genotype == "Control"] <- '#CC79A7'

map$DPW <- map$Age.Days - 21

otu <- read.csv("data/FD_OTU_mice21.csv",
                header = T, row.names = 1)

treats <-  c( "Mutant_Cohoused" , "Mutant_Separate",
              "Control_Cohoused", "Control_Separate")

#map.otu <- map.otu[ !is.na(map.otu$Phenotype.score) , ]

map <- map[!is.na(map$Phenotype.score), ]
map <- map[order(map$Age.Days) , ]

summary(table(map$Mouse.ID) > 1)

phenotype.lme <- lme4::lmer(map$Phenotype.score ~ map$Treatment.type + map$DPW + (1|map$Mouse.ID))
Treatment.null <- lme4::lmer(map$Phenotype.score ~ map$DPW + (1|map$Mouse.ID))
dpw.null <- lme4::lmer(map$Phenotype.score ~ map$Treatment.type+ (1|map$Mouse.ID))

anova(phenotype.lme, Treatment.null)
anova(phenotype.lme, dpw.null)

summary(phenotype.lme)



mice <- unique(map$Mouse.ID)

plot(map$Phenotype.score ~ map$DPW, pch = 21, bg = map$col, cex  =2,
     ylab = "Phenotype Score", xlab = "DPW", font = 3, font.lab = 4)

for(i in 1 : length(mice)){

  map.i <- map[map$Mouse.ID == mice[i] , ]

  if(nrow(map.i) > 1){

    print(i)

      for(j in 2: nrow(map.i)){

        segments(x0 =  map.i$DPW[j-1],
                 y0 =  map.i$Phenotype.score[j-1],
                 x1 =  map.i$DPW[j],
                 y1 =  map.i$Phenotype.score[j],
                 lwd = 2,
                 lty = 2,
                 col = map.i$col)
      }

  }
}

leg <- aggregate(map$col, list(map$discription), unique)
legend("bottomright", legend = leg$Group.1, pch = 21, pt.bg = leg$x)

########## compile source data file ###########33

Sup_data_6 <- cbind(map$ID.Gen.Tmt.Age,
                        map$Mouse.ID,
                        map$Treatment.type,
                        map$DPW,
                        map$Phenotype.score)

colnames(Sup_data_6) <- c("sample.id",
                          "mouse",
                          "treatment",
                          "DPW",
                          "Phenotype score")

Sup_data_6 <- cbind("", "", Sup_data_6, "", "")
Sup_data_6[1,1] <- "Supplemental Figure 6"

write.csv(Sup_data_6,"Sup_data_6.csv", row.names = F)


lme.res <- list(anova(phenotype.lme, Treatment.null),
                anova(phenotype.lme, dpw.null),
                summary(phenotype.lme))

names(lme.res) <- c("Effect of Treatment ANOVA", "Effect of DPW ANOVA", "Mixed Effect summary")
writeLines(capture.output(lme.res), "Sup_6_tests.txt")


