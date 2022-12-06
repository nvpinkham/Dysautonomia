
source("R/human_functions.22.11.27.R")

map <- read.csv("data/Map_human_microbiome.22.11.21.csv")
map.stool <- read.csv("data/Map_human_metabolome_stool.22.11.21.csv")
map.serum <- read.csv("data/Map_human_metabolome_serum.22.11.21.csv")


map$Sample.type <- "Stool 16S"
map.stool$Sample.type <- "Stool NMR"
map.serum$Sample.type <- "Serum NMR"

all <- merge(map, map.stool, all = T)
all <- aggregate(all, by = list(all$sample.id), unique.collapse)

all <- merge(all, map.serum, all = T)
all <- aggregate(all, by = list(all$sample.id), unique.collapse)

sup1a <- all[,2:4]

all$paired.relative <- gsub(" ", "", all$paired.relative)
pr <- strsplit(all$paired.relative, ",")
pr <- sapply(pr, unique)
sup1a$paired.relative <- sapply(pr, paste, collapse = ", ")
sup1a$paired.relative[sup1a$Disease.state != "Patient"] <- ""

sup1a$age <- all$age
sup1a$BMI <- all$BMI
sup1a$`Oral Foods (vs. exclusively eating with G tube)` <- all$Oral.foods#  ability to eat food orally (versus patients that exclusively ate using a gastrostomy (G-) tube
sup1a$`G tube (patient has gastrostomy but may still consume foods orally)` <- all$G.tube
sup1a$`Antibiotics (within previous 3 months)`<- all$Antibiotic

sup1a$ALP <- all$ALP
sup1a$ALT <- all$ALT
sup1a$AST <- all$AST
sup1a$BUN <- all$BUN
sup1a$Creatinine <- all$Creatinine
sup1a$eGFR <- all$eGFR

source("Figure_1/Figure_1_panel_A.22.11.29.R")
sup1a$`Included in Microbiome Ordination (Fig 1a)` <- all$sample.id %in% map$sample.id

source("Figure_1/Figure_1_Panel_E.22.11.28.R")
sup1a$`Included in Metabolome Stool Ordination (Fig 1e)` <- all$sample.id %in% map.stool$sample.id

source("Figure_1/Figure_1_panel_F.22.12.01.R")
sup1a$`Included in Metabolome Serum Ordination (Fig 1f)` <- all$sample.id %in% map.serum$sample.id

str(sup1a)

sup1a[sup1a == T] <- "yes"
sup1a[sup1a == F] <- "no"

write.csv(sup1a, "Sup_table_1a.csv")

################################################################################

pairs <- strsplit(all$pairs, ", ")
pairs <- unique(unlist(pairs))
pairs <- pairs[pairs != "NP"]

sup1b <- as.data.frame(pairs)


source("Figure_1/Figure_1_panel_C.22.11.29.R")
map.p <- map[map$Disease.state == "Patient", ]

pairs.filtered <- strsplit(map.p$pairs, ", ")
pairs.filtered <- unique(unlist(pairs.filtered))

cbind(pairs, pairs %in% pairs.filtered)

sup1b$`Pair Included in Alpha Diversity Analysis (Fig. 1 c, d)` <- pairs %in% pairs.filtered

source("Figure_2/Figure_2_Panel_B.22.11.28.R")
map.stool.p <- map.stool[map.stool$Disease.state == "Patient", ]
pairs.filtered1 <- strsplit(map.stool.p$pairs, ", ")
pairs.filtered1 <- unique(unlist(pairs.filtered1))

sup1b$`Included in Metabolite Stool Analysis (Fig. 2b)` <- pairs %in% pairs.filtered1

map.serum.p <- map.serum[map.serum$Disease.state == "Patient", ]

pairs.filtered2 <- strsplit(map.serum.p$pairs, ", ")
pairs.filtered2 <- unique(unlist(pairs.filtered2))

sup1b$`Included in Metabolite Serum Analysis (Fig. 2b)` <- pairs %in% pairs.filtered2


source("Figure_2/Figure_2_panel_C.22.11.29.R")
# Included in Targeted Stool TMA Analysis (Fig. 2c)
tma

tma.p <- tma[tma$Disease.state == "Patient", ]
pairs.filtered3 <- strsplit(tma.p$pairs, ", ")
pairs.filtered3 <- unique(unlist(pairs.filtered3))

cbind(pairs, pairs %in% pairs.filtered3)

sup1b$`Included in Targeted Stool TMA Analysis (Fig. 2c)` <- pairs %in% pairs.filtered3
#################
source("Figure_2/Figure_2_panel_D.22.11.29.R")
# Included in Targeted Stool TMA Analysis (Fig. 2c)
tmao

tmao.p <- tmao[tmao$Disease.state == "Patient", ]
pairs.filtered4 <- strsplit(tmao.p$pairs, ", ")
pairs.filtered4 <- unique(unlist(pairs.filtered4))

cbind(pairs, pairs %in% pairs.filtered4)

sup1b$`Included in Targeted Stool TMAO Analysis (Fig. 2d)` <- pairs %in% pairs.filtered4

sup1b[sup1b==T] <- "yes"
sup1b[sup1b==F] <- "no"
write.csv(sup1b, "Sup_table_1b.csv")





