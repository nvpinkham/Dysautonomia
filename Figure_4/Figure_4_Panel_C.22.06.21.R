
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

not.in.map <- otu[!row.names(otu) %in% map$ID.Gen.Tmt.Age,]
row.names(not.in.map) # these mice are not from the third generation
rm(not.in.map)

map.otu <- map[!is.na(map$X16SFile.ID) , ]
map.otu <- map.otu[map.otu$DOB > "2018-12-31" , ] # remove 2nd gen mice!
#  Mouse IDs 0-300 are third.


trans <- readxl::read_xlsx("data/GutTransit_MouseMasterData.xlsx")

trans$discription <- paste0(trans$Genotype, "_", trans$Treatment.type)

trans$discription <- factor(trans$discription,
                            levels = c("Control_Separate",
                                       "Control_Cohoused",
                                       "Mutant_Cohoused",
                                       "Mutant_Separate"))

sum(table(map.otu$Gut.Transit.Time.hours))

#setEPS()       # Set postscript arguments
#postscript(paste0("Figure4_panelC_", Sys.Date(), ".eps"),          width = 10, height = 6)

table(trans$discription)

vioplot(trans$Gut.Transit.Time.hours ~ trans$discription,
        drawRect = F,
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
t.test(trans$Gut.Transit.Time.hours[pick] ~ trans$discription[pick])

segments(2, 6, 3, 6)
text(2.5, 5.9, "n.s.", pos = 3, font = 4)

pick <- trans$Genotype == "Mutant"
t.res <- t.test(trans$Gut.Transit.Time.hours[pick] ~ trans$discription[pick])
t.res

segments(3, 6.6, 4, 6.6)
text(3.75, 6.75, paste("p =",
                      round(t.res$p.value, 4)), pos = 2, font = 4)


pick <- trans$discription %in% c("Control_Cohoused", "Mutant_Separate")
t.res <- t.test(trans$Gut.Transit.Time.hours[pick] ~ trans$discription[pick])
t.res

segments(2, 6.9, 4, 6.9)
text(3.25, 7.05, paste("p =",
                       round(t.res$p.value, 4)), pos = 2, font = 4)


########################

pick <- trans$discription %in% c("Control_Separate", "Mutant_Separate")
t.res <- t.test(trans$Gut.Transit.Time.hours[pick] ~ trans$discription[pick])
t.res

segments(1, 7.2, 4, 7.2)
text(2.75, 7.35, paste("p =",
                      round(t.res$p.value, 4)), pos = 2, font = 4)


#dev.off()

