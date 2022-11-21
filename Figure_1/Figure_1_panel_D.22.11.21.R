
# makes figure 1 panel C
# this should take around 15 minutes

source("R/human_functions.22.09.14.R")

map <- read.csv("data/Map_human_microbiome.22.11.21.csv")
row.names(map) <- map$sample.id

map <- read.csv("data/Map_human_microbiome.22.11.21.csv")
row.names(map) <- map$sample.id

dummy <- grep("[-]", map$Collection.date)
keep.dummy <- as.Date(map$Collection.date[dummy], format = "%Y-%m-%d")
map$Collection.date <- as.Date(map$Collection.date, format = "%m/%d/%y")# dummy collection dates to exclude from analysis
map$Collection.date[dummy] <- keep.dummy


otu <- read.csv("data/FD_OTU_human21.csv", row.names = 1)

reconcile(map, otu)
otu <- otu[!is.na(map$Collection.date) , ]
map <- map[!is.na(map$Collection.date) , ]

row.names(map) <- map$sample.id

filter <- NULL
for(i in 1:length(map$sample.id)){
  # remove samples that can't be paired

  map.i <- map[i,]
  fam.i <- map[map$Family.ID == map.i$Family.ID,]

  pairs.pos <- fam.i[fam.i$Disease.state != map.i$Disease.state , ] #possible pairs

  if(nrow(pairs.pos) > 0){

    closest.date <- min(abs(pairs.pos$Collection.date - map.i$Collection.date))
    filter[i] <- closest.date < 90
  }else{
    filter[i] <- FALSE

  }
}

otu <- otu[filter,]
map <- map[filter,]

reconcile(map, otu)


map$prop.shared <- NA

prop.shared <- matrix(nrow = nrow(map), ncol = 999)
row.names(prop.shared) <- map$sample.id
prop.unique <- prop.shared


t.res.shared.stat <- t.res.shared.p <- t.res.shared.df <- t.res.shared.ef <- NULL

p.shared <- list() # only in one group
p.shared <- p.shared[1:999]
r.shared <- p.shared
pats <- unique(map$DonorID[map$Disease.state == "Patient"])


#proportion of population
for(i in  1 : 999){

  print(i)

  map.i <- paired.map(map, pats)
  otu.i <- otu[match(map.i$sample.id, row.names(otu)) , ]

  for(j in 1 : length(pats)){

    fam.j <- map.i$Family.ID[map.i$DonorID == pats[j]]

    p <- otu.i[map.i$Disease.state == "Patient" & map.i$Family.ID == fam.j, ]
    r <- otu.i[map.i$Disease.state != "Patient" & map.i$Family.ID == fam.j, ]

    p <- p[ , colSums(p) > 0]
    r <- r[ , colSums(r) > 0]

    p.shared[[i]][j] <- sum(p[ colnames(p) %in% colnames(r)]) / 100# percent of microbiome (rared to 10000)
    r.shared[[i]][j] <- sum(r[ colnames(r) %in% colnames(p)]) / 100
  }

  prop.shared[match(row.names(map.i), row.names(prop.shared)) , i] <- c(p.shared[[i]],
                                                                        r.shared[[i]])

  t.res.shared.res <- t.test(p.shared[[i]], r.shared[[i]], paired = T)
  t.res.shared.p[i] <-  t.res.shared.res$p.value
  t.res.shared.ef[i] <- t.res.shared.res$estimate
  t.res.shared.df[i] <- t.res.shared.res$parameter
  t.res.shared.stat[i] <- t.res.shared.res$statistic



}

map$prop.shared <- rowMeans(prop.shared, na.rm = T)

saveRDS(t.res.shared.p, "Microbiome_Permution_Analysis/t.res.shared.p.rds")
saveRDS(t.res.shared.df, "Microbiome_Permution_Analysis/t.res.shared.df.rds")
saveRDS(t.res.shared.ef, "Microbiome_Permution_Analysis/t.res.shared.ef.rds")
saveRDS(t.res.shared.stat, "Microbiome_Permution_Analysis/t.res.shared.stat.rds")

t.res.shared.p <- readRDS("Microbiome_Permution_Analysis/t.res.shared.p.rds")
t.res.shared.df <- readRDS("Microbiome_Permution_Analysis/t.res.shared.df.rds")
t.res.shared.ef <- readRDS("Microbiome_Permution_Analysis/t.res.shared.ef.rds")

set.seed(42)
map$dummy <- 1
map$dummy[map$Disease.state != "Patient" ] <- 2
map$dummy  <- map$dummy + sample(-100 : 100, nrow(map), replace = T)/2000

map.pat <- map[map$Disease.state == "Patient", ]
map.rel <- map[map$Disease.state == "Relative", ]


################################################################################
#####                             SHARED                         ###############


vioplot(map.pat$prop.shared,
        map.rel$prop.shared,
        drawRect = F, col=unique(map$col),
        names = c("Patients", "Relatives"),
        ylab = "% shared OTUs")

segments(.7, mean(map.pat$prop.shared),
         1.3, mean(map.pat$prop.shared), lwd = 2)

segments(1.7, mean(map.rel$prop.shared),
         2.3, mean(map.rel$prop.shared), lwd = 2)


families <- unique(map$Family.ID)

for(i in 1: length(families)){

  map.pick <- map.pat[map.pat$Family.ID ==  families[i] , ]
  map.pick2 <- map.rel[map.pat$Family.ID ==  families[i] , ]

  segments(y0 = map.pick$prop.shared, x0 = map.pick$dummy,
           y1 = map.pick2$prop.shared, x1 = map.pick2$dummy,
           col = 8)

  points(map.pick$dummy, map.pick$prop.shared, pch = 21, bg = 8)
  points(map.pick2$dummy, map.pick2$prop.shared, pch = 21, bg = 8)

}

title(paste("Paired T test mean p value =",  round(mean(t.res.shared.p, na.rm = T), 5), "\ndf = ",
            paste(unique(t.res.shared.df), collapse = " "),
            "\nmean of the differences", round(mean(t.res.shared.ef, na.rm = T), 5)))

################################################################################

source_data.1d <- cbind(map$sample.id,  map$DonorID, map$Family.ID, map$Disease.state, map$prop.shared)
colnames(source_data.1d) <- c("sample.id",  "DonorID", "Family.ID", "Disease.state", "Mean percent shared with paired sample")

source_data.1d <- cbind("", source_data.1d)
source_data.1d[1,1] <- "Figure 1d"

write.csv(source_data.1d,"source_data_1d.csv")



