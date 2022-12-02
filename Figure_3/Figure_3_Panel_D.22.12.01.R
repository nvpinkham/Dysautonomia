

source("R/mice_functions.22.12.01.R")


dir.create("processing")

meta.6 <-  t(read.csv("data/6mo2018mousestool_umolepermg.csv",
                      header = T, row.names = 1))
meta.6 <- as.data.frame(meta.6)

meta.6 <- norm(row.names(meta.6), meta.6)

#########

gtreatment_group <- meta.6$Label

meta.6$Label <- NULL

map.6 <- as.data.frame(cbind(row.names(meta.6),
                             gtreatment_group,
                             as.numeric(as.factor(gtreatment_group))))
colnames(map.6) <- c("id", "discription", "col")

map.6$col[map.6$discription == "FD_6mo"] <-  "#D55E00"
map.6$col[map.6$discription != "FD_6mo"] <-  "#F0E442"

plot.meta(meta.pick = meta.6, map.pick = map.6)

########## compile source data file ###########33
all(row.names(meta.6) == map.6$id)
map.6$genotype <- gsub("_6mo", "", map.6$discription)

source_data_3d <- cbind(map.6$id, map.6$genotype)
colnames(source_data_3d) <- c("sample.id", "genotype")
source_data_3d <- cbind("", "", source_data_3d, "", "")
source_data_3d[1,1] <- "figure 3d"

source_data_3d[1,6] <- "Metabolite concentrations"

source_data_3d <- cbind(source_data_3d, as.matrix(meta.6))

write.csv(source_data_3d,"source_data_3d.csv", row.names = F)
set.seed(42)
res <- adonis2(dist(meta.6) ~ gtreatment_group, permutations = 9999)
res

res <- capture.output(res)
writeLines(res, "Statistical_summaries/Fig_3d_tests.txt")

