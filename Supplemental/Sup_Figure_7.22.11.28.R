
source("R/mice_functions.22.11.27.R")

trans <- readxl::read_xlsx("data/GutTransit_MouseMasterData.xlsx")

trans$DPW <- trans$Age.Days -21

trans$discription <- paste0(trans$Genotype, "_", trans$Treatment.type)

trans$discription <- factor(trans$discription,
                            levels = c("Control_Separate",
                                       "Control_Cohoused",
                                       "Mutant_Cohoused",
                                       "Mutant_Separate"))

trans$col <- NA
trans$col[trans$Treatment.type =="Separate" & trans$Genotype == "Mutant" ] <- '#D55E00'
trans$col[trans$Treatment.type =="Separate" & trans$Genotype == "Control"] <- '#F0E442'
trans$col[trans$Treatment.type =="Cohoused" & trans$Genotype == "Mutant"] <-  '#009E73'
trans$col[trans$Treatment.type =="Cohoused" & trans$Genotype == "Control"] <- '#CC79A7'



plot(trans$Gut.Transit.Time.hours ~ trans$DPW,
     cex = 3,
     pch =21,
     bg = trans$col,
     xlab = "DPW",
     ylab = "Gut Transit Time (hours)")

leg <- aggregate(trans$col, list(trans$discription), unique)
legend("bottomright", legend = leg$Group.1, pch = 21, pt.bg = leg$x)

treats <- unique(trans$discription)

for(i in 1:4){

  pick <- trans$discription == treats[i]

  res <- cor.test(trans$Age.Days[pick], trans$Gut.Transit.Time.hours[pick], method = "spearman")

  print(res)
}


########################

pick <- trans$discription %in% c("Control_Separate", "Mutant_Separate")
t.res <- t.test(trans$Gut.Transit.Time.hours[pick] ~ trans$discription[pick])
t.res


library(lme4)

a <- sapply(strsplit(trans$ID.Gen.Tmt.Age, "[.]"), "[[", 1)
sort(table(a))# only done once

trans <- trans[trans$Genotype == "Mutant",]
trans$Cage.ID[is.na(trans$Cage.ID)] <- "unknown"

Transit.lme <- lmer(trans$Gut.Transit.Time.hours ~ trans$discription + trans$DPW + (1|trans$Cage.ID),
                    data = trans)
DPW.null <- lmer(trans$Gut.Transit.Time.hours ~ trans$discription + (1|trans$Cage.ID),
                 data = trans)

anova(Transit.lme, DPW.null)

row.names(trans)
Sup_data_7 <- cbind(trans$ID.Gen.Tmt.Age,
                        trans$Treatment.type,
                        trans$Age.Days - 21,
                        trans$Gut.Transit.Time.hours)

colnames(Sup_data_7) <- c("sample.id",
                              "treatment",
                              "DPW",
                              "Gut Transit Time (hours)")

Sup_data_7 <- cbind("", "", Sup_data_7)
Sup_data_7[1,1] <- "Supplemental figure 7"

write.csv(Sup_data_7,"Sup_data_7.csv", row.names = F)

res <- list(res, summary(Transit.lme))
names(res) <- c("DPW ANOVA", "model Summary")
writeLines(as.character(res), "Sup_7_test.txt")



