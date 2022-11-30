
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


gt.full <- lm(trans$Gut.Transit.Time.hours ~ trans$Genotype +  trans$Treatment.type + trans$DPW)

lm.sum <- summary.lm(gt.full)
lm.sum
confint(gt.full) # DPW not significant

trans2 <- trans[trans$Genotype == "Mutant",]
gt.mu <- lm(trans2$Gut.Transit.Time.hours ~ trans2$Treatment.type + trans2$DPW)

summary.lm(gt.mu)
confint(gt.mu) #


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

write.csv(Sup_data_7,"source_data/Sup_data_7.csv", row.names = F)

res <- list(lm.sum,
            confint(gt.full))


names(res) <- c("LM summary", "CI")

writeLines(capture.output(res),
           "Statistical_summaries/Sup_7_tests.txt")

