

source("R/mice_functions.22.11.27.R")

map <- read.csv("data/Venn_map.22.11.25.csv", header = T, row.names = 1)
otu <- read.csv("data/Venn_otu.22.11.25.csv", header = T, row.names = 1)

labdsv::reconcile(map, otu)
otu <- vegan::rrarefy(otu, 5000)

venn2otu(otu, map, "mouse gut metagenome", "human gut metagenome", min.sites = 5)

############################ Write Source data #################################

sample.id <- map$Experiment.Title
Organism <- map$Organism.Name

sup_data_9 <- cbind("", "", sample.id, Organism, "", "",  as.matrix(otu))
sup_data_9[1,1] <- "Supplemental figure 9"
sup_data_9[1,2] <- "Combined human and mice analysis"
sup_data_9[1,6] <- "OTU counts"

write.csv(sup_data_9,"sup_data_9.csv", row.names = F)
################################################################################




