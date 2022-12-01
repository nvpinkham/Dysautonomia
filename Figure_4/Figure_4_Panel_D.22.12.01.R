
source("R/mice_functions.22.12.01.R")

trans <- readxl::read_xlsx("data/GutTransit_MouseMasterData.xlsx")

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

#setEPS()       # Set postscript arguments
#postscript(paste0("Figure4_panelC_", Sys.Date(), ".eps"),          width = 10, height = 6)

table(trans$discription)

vioplot::vioplot(trans$Gut.Transit.Time.hours ~ trans$discription,
                 drawRect = F,
                 #col = c("#56B4E9", "#CC79A7", "#009E73", "#E69F00"),
                 col = c("#F0E442", "#CC79A7", "#009E73", "#D55E00"),
                 xlim = c(.5, 5.5),
                 ylim = c(2, 7.5),

                 names = paste0(levels(trans$discription),
                                "\nn=",
                                table(trans$discription)),

                 xlab = "",
                 ylab = "Gut Transit Time (hours)")

agg <- aggregate(trans$Gut.Transit.Time.hours, list( trans$discription), mean)
agg$Group.1

segments(.7, agg$x[3], 1.3,  agg$x[3], lwd = 2)
segments(1.7, agg$x[2], 2.3,  agg$x[2], lwd = 2)
segments(2.7, agg$x[1], 3.3,  agg$x[1], lwd = 2)
segments(3.7, agg$x[4], 4.3,  agg$x[4], lwd = 2)


trans$Age.Days.norm <- scale( trans$Age.Days)
trans$Age.Days.norm <- trans$Age.Days.norm + diff(c(min(trans$Age.Days.norm, na.rm = T), 1.25))

#trans$Age.Days.norm <- trans$Age.Days.norm * 2

trans <- trans[rev(order(trans$Age.Days)), ]

points(trans$Gut.Transit.Time.hours ~ trans$discription,
       pch = 21, bg = 8, cex = trans$Age.Days.norm)

a <- trans$Age.Days.norm
b <- trans$Age.Days

ages <- trans$Age.Days[!is.na(trans$Gut.Transit.Time.hours )]

s.lm <- lm(a ~ b)


#nd <- data.frame(b = c(100, 200, 300))
nd <- data.frame(b = c(min(ages), 300, max(ages)))

predict(s.lm , newdata = nd)

legend(pch = 21,
       pt.bg  = 8,
       cex = 1,
       bty = "n",
       pt.cex = predict(s.lm , newdata = nd)  ,
       legend = paste("\n ", nd$b, "\n") ,
       5, 7)

text(5.15, 6.85, "Age\n(Days)")

rect(1.5, 2 , 3.5, 5.55, lty = 2)

pick <- trans$Treatment.type == "Cohoused"
cohoused.t <- t.test(trans$Gut.Transit.Time.hours[pick] ~ trans$discription[pick])

segments(2, 6, 3, 6)
text(2.5, 5.9, "n.s.", pos = 3, font = 4)

pick <- trans$Genotype == "Mutant"
mutant.t <- t.test(trans$Gut.Transit.Time.hours[pick] ~ trans$discription[pick])
mutant.t

segments(3, 6.6, 4, 6.6)
text(3.75, 6.75, paste("p =",
                      round(mutant.t$p.value, 4)), pos = 2, font = 4)


pick <- trans$discription %in% c("Control_Cohoused", "Mutant_Separate")
Control_Cohoused_Mutant_Separate.t <- t.test(trans$Gut.Transit.Time.hours[pick] ~ trans$discription[pick])
Control_Cohoused_Mutant_Separate.t

segments(2, 6.9, 4, 6.9)
text(3.25, 7.05, paste("p =",
                       round(Control_Cohoused_Mutant_Separate.t$p.value, 4)),
     pos = 2, font = 4)

########## compile source data file ###########33

row.names(trans)
source_data_4d <- cbind(trans$ID.Gen.Tmt.Age,
                        trans$Treatment.type,
                        trans$Age.Days - 21,
                        trans$Gut.Transit.Time.hours)

colnames(source_data_4d) <- c("sample.id",
                              "treatment",
                              "DPW",
                              "Gut Transit Time (hours)")

source_data_4d <- cbind("", "", source_data_4d)
source_data_4d[1,1] <- "figure 4d"

write.csv(source_data_4d,"source_data_4d.csv", row.names = F)

#################
unique(trans$discription)

pick <- trans$Genotype == "Control"
control.t <- t.test(trans$Gut.Transit.Time.hours[pick] ~ trans$discription[pick])

pick <- trans$discription %in% c("Control_Separate", "Mutant_Cohoused")
Control_Separate_Mutant_Cohouse.t <- t.test(trans$Gut.Transit.Time.hours[pick] ~ trans$discription[pick])

pick <- trans$discription %in% c("Control_Separate", "Mutant_Separate")
Control_Separate_Mutant_Separate.t <- t.test(trans$Gut.Transit.Time.hours[pick] ~ trans$discription[pick])

t.res <- list(cohoused.t,
              mutant.t,
              Control_Cohoused_Mutant_Separate.t,
              control.t,
              Control_Separate_Mutant_Cohouse.t,
              Control_Separate_Mutant_Separate.t)

names(t.res) <- c("Cohoused mutants and controls",
                  "Mutants",
                  "Control_Cohoused and Mutant_Separate",
                  "Controls",
                  "Control_Separate and Mutant_Cohouse",
                  "Control_Separate and Mutant_Separate")

ps <- c(t.res[[1]]$p.value,
        t.res[[2]]$p.value,
        t.res[[3]]$p.value,
        t.res[[4]]$p.value,
        t.res[[5]]$p.value,
        t.res[[6]]$p.value)

ps.fdr <- p.adjust(ps)

for(i in 1 : 6){
  t.res[[i]]$method <- paste(t.res[[i]]$method, "; p after FDR = ", ps.fdr[i])
}

t.res[[7]] <- table(trans$discription, trans$Sex)

writeLines(capture.output(t.res), "Statistical_summaries/Fig_4d_tests.txt")


