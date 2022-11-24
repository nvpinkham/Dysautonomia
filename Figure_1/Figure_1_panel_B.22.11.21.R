source("R/human_functions.22.11.21.R")

# import data
map <- read.csv("data/Map_human_microbiome.22.11.21.csv")
row.names(map) <- map$sample.id
otu <- read.csv("data/FD_OTU_human21.csv", row.names = 1)

labdsv::reconcile(map, otu) # make sure samples match

otu.bc <- vegdist(otu) # calc dissimilarity

between <- as.matrix(otu.bc)[map$Disease.state == "Patient",
                             map$Disease.state != "Patient"]

between <- as.vector(between)

otu.pat <- otu[map$Disease.state == "Patient",]
pats <- vegdist(otu.pat)

otu.rel <- otu[map$Disease.state != "Patient",]
rels <- vegdist(otu.rel)


vioplot(pats, between, rels,
        ylab = "BC distances",
        col = c("#E69F00", "grey","#56B4E9"),
        names = c("within Patients",
                  "between Patients and Relatives",
                  "within Relatives"))


t.test(pats, rels)
t.test(pats, between)

pats <- as.vector(pats)
between <- as.vector(between)
rels <- as.vector(rels)

length(between)

pats <- c(pats, rep("", length(between) - length(pats)))
rels <- c(rels, rep("", length(between) - length(rels)))

source_data.1b <- cbind(pats, between, rels)
colnames(source_data.1b) <- c("BC distance within patients",
                              "BC distance between patients and relatives",
                              "BC distance within Relatives")

source_data.1b <- cbind("", source_data.1b)
source_data.1b[1,1] <- "Figure 1b"

write.csv(source_data.1b,"source_data_1b.csv")

