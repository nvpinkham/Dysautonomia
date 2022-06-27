
source("R/human_functions.22.06.18.R")

meta.serum <- read.csv("data/HumanSerumUnpaired_biomassnorm_metalabels.csv", row.names = 1)
map.serum <- read.csv("data/Map_human_metabolome_serum.csv", row.names = 1)

meta.stool <- read.csv("data/HumanStoolUnpaired_biomassnorm_metalabels.csv", row.names = 1)
map.stool <- read.csv("data/Map_human_metabolome_stool.csv", row.names = 1)

map.serum$Collection.date <- as.Date(map.serum$Collection.date )
map.stool$Collection.date <- as.Date(map.stool$Collection.date )

tmao <- map.serum[ !is.na(map.serum$tmao) , ]
table(tmao$Family.ID)

tma <- map.stool[ !is.na(map.stool$tma) , ]
table(tma$Family.ID)

tmao$tmao.ln <- log(tmao$tmao)
tma$tma.ln <- log(tma$tma)

tma <- tma[order(tma$Family.ID), ]
tma$Family.ID

table(tma$Disease.state)
table(tma$Family.ID)

tma.36 <- tma[tma$Family.ID == 36 , ]
tma.36$Collection.date

# only pair with one sample
tma$Family.ID[tma$sample.id %in% c("FD082", "FD083")] <- 360
tmao$Family.ID[tmao$sample.id %in% c("FD082", "FD083")] <- 360

pats <- tma$DonorID[tma$Disease.state == "Patient"]
tma2 <- paired.map(tma, pats)

qqnorm(tma2$tma[tma2$Disease.state == "Patient"])
qqnorm(tma2$tma[tma2$Disease.state != "Patient"])

shapiro.test(tma2$tma[tma2$Disease.state == "Patient"])
shapiro.test(tma2$tma[tma2$Disease.state != "Patient"])

wilcox.test(tma2$tma ~ tma2$Disease.state, paired = T)
t.test(tma2$tma.ln ~ tma2$Disease.state, paired = T)

qqnorm(tma2$tma.ln[tma2$Disease.state == "Patient"])
qqnorm(tma2$tma.ln[tma2$Disease.state != "Patient"])

shapiro.test(tma2$tma.ln[tma2$Disease.state == "Patient"])
shapiro.test(tma2$tma.ln[tma2$Disease.state != "Patient"])

t.res <- t.test(tma2$tma.ln ~ tma2$Disease.state, paired= T)

paired.violin.beta(map = tma2, var = tma2$tma.ln, plot.each = F)
mtext(side = 2,at =  2, line = 2, "Natural log uM", font = 4)
title(paste0("t-test\np-value =", round(t.res$p.value, 5),
             ", effect size= ", round(t.res$estimate, 5)),
      font = 4)

t.test(tmao$tmao ~ tmao$Disease.state, paired= T)
t.test(tmao$tmao.ln ~ tmao$Disease.state, paired= T)


all(tmao$sample.id %in% map.serum$sample.id)
all(tma$sample.id %in% map.stool$sample.id)


map.tmao$tmao.ln <- log(tmao$tmao)

tmao.pick$sample.id
tmao$sample.id



t.res <- t.test(tma2$tma.ln ~ tma2$Disease.state, paired= T)

paired.violin.beta(map = tma2, var = tma2$tma.ln, plot.each = F)
mtext(side = 2,at =  2, line = 2, "Natural log uM", font = 4)
title(paste0("t-test\np-value =", round(t.res$p.value, 5),
             ", effect size= ", round(t.res$estimate, 5)),
      font = 4)


t.res <- t.test(tmao$tmao.ln ~ tmao$Disease.state, paired= T)

paired.violin.beta(map = tmao, var = tmao$tmao.ln, plot.each = F)
mtext(side = 2,at =  0, line = 2, "Natural log uM", font = 4)
title(paste0("t-test\np-value =", round(t.res$p.value, 5),
             ", effect size= ", round(t.res$estimate, 5)),
      font = 4)

