
source("R/mice_functions.22.06.22.R")

map <- read.csv("data/Map_CoHoused_mice.8.16.csv")
map$discription <- paste0(map$Genotype, "_", map$Treatment.type)
map$ID.Gen.Tmt.Age <- paste0("X", map$ID.Gen.Tmt.Age)
map$DOB <- as.Date(map$DOB , "%m/%d/%y")

map$X16SFile.ID[map$X16SFile.ID == "" ] <- NA

map$col[map$Treatment.type =="Separate" & map$Genotype == "Mutant" ] <- '#D55E00'
map$col[map$Treatment.type =="Separate" & map$Genotype == "Control"] <- '#F0E442'
map$col[map$Treatment.type =="Cohoused" & map$Genotype == "Mutant"] <-  '#009E73'
map$col[map$Treatment.type =="Cohoused" & map$Genotype == "Control"] <- '#CC79A7'

otu <- read.csv("data/FD_OTU_mice21.csv",
                header = T, row.names = 1)

treats <-  c( "Mutant_Cohoused" , "Mutant_Separate",
              "Control_Cohoused", "Control_Separate")

map.1 <- get.close(map, 100)
t1 <- t.test(map.1$Phenotype.score ~ map.1$Treatment.type)

map.2 <- get.close(map, 200)
t2 <- t.test(map.2$Phenotype.score ~ map.2$Treatment.type)

map.3 <- get.close(map, 300)
t3 <- t.test(map.3$Phenotype.score ~ map.3$Treatment.type)


map <- rbind(map.1,
             map.2,
             map.3)



map <- map[order(map$Age.Days) ,  ]

map$test <- paste("Day", map$Age.Bin, "\n", map$Treatment.type)
map$test <- as.factor(map$test)


#setEPS()       # Set postscript arguments
#postscript(paste0("Figure4_panelD_", Sys.Date(), ".eps"),
#           width = 10, height = 6)

vioplot::vioplot(map$Phenotype.score ~ map$test,
                 ylab = "Phenotype Score",
                 xlab = "",
                 ylim = c(0, 11),
                 drawRect = F,
                 # names = c("3 months\nCohoused", "3 months\nSeperated",
                 #          "6 months\nCohoused", "6 months\nSeperated",
                 #         "9 months\nCohoused", "9 months\nSeperated"),
                 col = c("#009E73", "#D55E00"))

mice <- sort(unique(map$Mouse.ID))


map$dummy <- as.numeric( map$test) + rnorm(length(map$test), mean = 0, sd = .03)

table(map$Mouse.ID, is.na(map$Phenotype.score))

points(map$Phenotype.score + rnorm(length(map$test), mean = 0, sd = .03) ~ map$dummy, col = map$col, pch = 21, bg = 8)


segments(1, 8.5, 2, 8.5)
text(1.5, 8.7, paste("p = ", round(t1$p.value, 4)), font = 4)
segments(3, 9.5, 4, 9.5)
text(3.5, 9.7, paste("p = ", round(t2$p.value, 4)), font = 4)
segments(5, 10.5, 6, 10.5)
text(5.5, 10.7, paste("p = ", round(t3$p.value, 4)), font = 4)

mean(map$Phenotype.score, na.rm = T)

agg <- aggregate(map$Phenotype.score, list( map$test), function(x) mean(x, na.rm=TRUE))

segments(.7, agg$x[1], 1.3, agg$x[1], lwd = 2)
segments(1.7, agg$x[2], 2.3, agg$x[2], lwd = 2)
segments(2.7, agg$x[3], 3.3, agg$x[3], lwd = 2)
segments(3.7, agg$x[4], 4.3, agg$x[4], lwd = 2)
segments(4.7, agg$x[5], 5.3, agg$x[5], lwd = 2)
segments(5.7, agg$x[6], 6.3, agg$x[6], lwd = 2)


#dev.off()

table(map$Age.Bin[!is.na(map$Phenotype.score)],
      map$discription[!is.na(map$Phenotype.score)])


########## compile source data file ###########33
source_data_4c <- cbind(map$ID.Gen.Tmt.Age,
                        map$Mouse.ID,
                        map$Treatment.type,
                        map$Age.Days - 21,
                        map$Phenotype.score)

colnames(source_data_4c) <- c("sample.id",
                              "mouse",
                              "treatment",
                              "DPW",
                              "Phenotype.score")

source_data_4c <- cbind("", "", source_data_4c)
source_data_4c[1,1] <- "figure 4c"

write.csv(source_data_4c,"source_data_4c.csv", row.names = F)












