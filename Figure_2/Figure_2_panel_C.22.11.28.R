
source("R/human_functions.22.11.27.R")

map.stool <- read.csv("data/Map_human_metabolome_stool.22.11.21.csv", row.names = 1)
map.stool$Collection.date <- as.Date(map.stool$Collection.date, format = "%m/%d/%y")# excel messes up dates by careful
map.stool$sample.id <- row.names(map.stool)

tma <- map.stool[ !is.na(map.stool$tma) , ]
table(tma$Family.ID)

tma$tma.ln <- log(tma$tma)

tma <- tma[order(tma$Family.ID), ]

table(tma$Disease.state)
table(tma$Family.ID)

tma.36 <- tma[tma$Family.ID == 36 , ]

# only pair with one sample
tma$Family.ID[tma$sample.id %in% c("FD082", "FD083")] <- 360

pats <- tma$DonorID[tma$Disease.state == "Patient"]
tma2 <- paired.map(tma, pats)

#setEPS()       # Set postscript arguments
#postscript(paste0("figure_2_panel_C_", Sys.Date(), ".eps"), width = 4, height = 3)

t.res <- t.test(tma2$tma.ln ~ tma2$Disease.state, paired= T)

paired.violin.beta(map = tma2, var = tma2$tma.ln, plot.each = F)
mtext(side = 2,at =  2, line = 2, "Natural log uM", font = 4)

title.t <- paste("TMA t-test\neffect size =",
                 round(diff(t.res$estimate), 5),
                 "p val = ", round(t.res$p.value, 5),
                 "\n95% CI =", round(t.res$conf.int[1], 5),
                 ",", round(t.res$conf.int[2], 5),
                 "\ndf = ", round(t.res$parameter, 5))

title(title.t,
      font = 4)

#dev.off()

source_data_2c <- cbind(tma2$sample.id,
                        tma2$Disease.state,
                        tma2$tma.ln)

colnames(source_data_2c) <- c("sample.id",
                              "Disease.state",
                              "natural log tma")

source_data_2c <- cbind("", "", source_data_2c)

source_data_2c[1,2] <- "Figure 2c"

write.csv(source_data_2c,"source_data/source_data_2c.csv", row.names = F)
t.res

t.res <- capture.output(t.res)
writeLines(t.res, "Statistical_summaries/Fig_2c_tests.txt")

