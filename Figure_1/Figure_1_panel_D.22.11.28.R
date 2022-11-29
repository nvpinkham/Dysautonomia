
# this should take around 15 minutes

source("R/human_functions.22.11.27.R")

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


#t.res.shared.stat <- t.res.shared.p <- t.res.shared.df <- t.res.shared.ef <- NULL

shared.t <- as.data.frame(matrix(nrow = 999, ncol = 6))
colnames(shared.t) <- c("P_val", "df", "effect_size", "t.stat", "conf_int_lo",  "conf_int_hi")


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

  res <- t.test(p.shared[[i]], r.shared[[i]], paired = T)

  shared.t$P_val[i] <- res$p.value
  shared.t$df[i] <- res$parameter
  shared.t$effect_size[i] <- res$estimate
  shared.t$t.stat[i] <- res$statistic
  shared.t$conf_int_lo[i] <- res$conf.int[1]
  shared.t$conf_int_hi[i] <- res$conf.int[2]


}

map$prop.shared <- rowMeans(prop.shared, na.rm = T)

saveRDS(shared.t, "Microbiome_Permution_Analysis/shared.t.rds")

shared.t <- readRDS("Microbiome_Permution_Analysis/shared.t.rds")


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

title(paste("Paired T test mean p value =",  round(mean(shared.t$P_val, na.rm = T), 5),
            "\nmean effect size =", round(mean(shared.t$effect_size, na.rm = T), 5),
            "\ndf =", round(mean(shared.t$df, na.rm = T), 5),
            "\nmean 95% CI =", round(mean(shared.t$conf_int_lo, na.rm = T), 2),"to",
            round(mean(shared.t$conf_int_hi, na.rm = T), 2)))
################################################################################

source_data.1d <- cbind(map$sample.id,  map$DonorID, map$Family.ID, map$Disease.state, map$prop.shared)
colnames(source_data.1d) <- c("sample.id",  "DonorID", "Family.ID", "Disease.state", "Mean percent shared with paired sample")

source_data.1d <- cbind("", source_data.1d)
source_data.1d[1,1] <- "Figure 1d"

write.csv(source_data.1d,"source_data_1d.csv")

shared.t.means <- apply(shared.t, 2, mean)
names(shared.t.means) <- paste("mean over 999 permutations", names(shared.t.means))
shared.t.means <- capture.output(shared.t.means)
writeLines(shared.t.means, "Fig_1d_tests.txt")


