
source("R/human_functions.22.11.27.R")

map.serum <- read.csv("data/Map_human_metabolome_serum.csv", row.names = 1)
map.serum$Collection.date <- as.Date(map.serum$Collection.date, format = "%Y-%m-%d")

tmao <- map.serum[ !is.na(map.serum$tmao) , ]
table(tmao$Family.ID)

tmao$tmao.ln <- log(tmao$tmao)

# only pair with one sample
tmao$Family.ID[tmao$sample.id %in% c("FD082", "FD083")] <- 360


#setEPS()       # Set postscript arguments
#postscript(paste0("figure_2_panel_D_", Sys.Date(), ".eps"), width = 4, height = 3)

t.res <- t.test(tmao$tmao.ln ~ tmao$Disease.state, paired= T)

paired.violin.beta(map = tmao, var = tmao$tmao.ln, plot.each = F)
mtext(side = 2,at =  0, line = 2, "Natural log uM", font = 4)


title.t <- paste("TMAO t-test\neffect size =",
                 round(diff(t.res$estimate), 5),
                 "p val = ", round(t.res$p.value, 5),
                 "\n95% CI =", round(t.res$conf.int[1], 5),
                 ",", round(t.res$conf.int[2], 5),
                 "\ndf = ", round(t.res$parameter, 5))

title(title.t,
      font = 4)

#dev.off()

source_data_2d <- cbind(tmao$sample.id,
                        tmao$Disease.state,
                        tmao$tmao.ln)

colnames(source_data_2d) <- c("sample.id",
                              "Disease.state",
                              "natural log tmao")

source_data_2d <- cbind("", "", source_data_2d)

source_data_2d[1,2] <- "Figure 2d"

write.csv(source_data_2d,"source_data_2d.csv", row.names = F)
t.res
