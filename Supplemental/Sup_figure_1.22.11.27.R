

source("R/human_functions.22.11.27.R")
# run time : 10-20 minutes
num.perm = 9

map <- read.csv("data/Map_human_microbiome.22.11.21.csv")
row.names(map) <- map$sample.id
dummy <- grep("[-]", map$Collection.date)
keep.dummy <- as.Date(map$Collection.date[dummy], format = "%Y-%m-%d")
map$Collection.date <- as.Date(map$Collection.date, format = "%m/%d/%y")# dummy collection dates to exclude from analysis
map$Collection.date[dummy] <- keep.dummy

otu <- read.csv("data/FD_OTU_human21.csv", row.names = 1)

table(map$Disease.state)
length(unique(map$Family.ID))

reconcile(map, otu)

map$Collection.date <- as.Date(map$Collection.date)
otu <- otu[!is.na(map$Collection.date) , ]
map <- map[!is.na(map$Collection.date) , ]

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

row.names(map) <- map$sample.id
reconcile(map, otu)


prop.inv <- matrix(nrow = nrow(map), ncol = num.perm)
row.names(prop.inv) <- map$sample.id
prop.unique <- prop.inv

num.perm = 999

t.inv.lo <- t.inv.hi <- t.inv.stat <- t.inv.p <- t.inv.df <- t.inv.ef <- NULL

pats <- unique(map$DonorID[map$Disease.state == "Patient"])

cor.inv.res <- NULL
for(i in  1 : num.perm){
  # paired permutation analysis
  # Each iteration of loop is a permutation

  print(i)

  map.i <- paired.map(map, pats) # randomly sampled each time

  t.inv.res <- t.test( map.i$invsimp[map.i$Disease.state == "Patient"],
                       map.i$invsimp[map.i$Disease.state == "Relative"],
                       paired = T)

  t.inv.p[i] <-  t.inv.res$p.value
  t.inv.ef[i] <- t.inv.res$estimate
  t.inv.df[i] <- t.inv.res$parameter
  t.inv.stat[i] <- t.inv.res$statistic
  t.inv.lo[i] <- t.inv.res$conf.int[1]
  t.inv.hi[i] <- t.inv.res$conf.int[2]
}

# save results from permutation analysis
dir.create("Microbiome_Permution_Analysis")

saveRDS(t.inv.p, "Microbiome_Permution_Analysis/t.inv.p.rds")
saveRDS(t.inv.df, "Microbiome_Permution_Analysis/t.inv.df.rds")
saveRDS(t.inv.ef, "Microbiome_Permution_Analysis/t.inv.ef.rds")
saveRDS(t.inv.stat, "Microbiome_Permution_Analysis/t.inv.stat.rds")


t.inv.p <- readRDS("Microbiome_Permution_Analysis/t.inv.p.rds")
t.inv.df <- readRDS("Microbiome_Permution_Analysis/t.inv.df.rds")
t.inv.ef <- readRDS("Microbiome_Permution_Analysis/t.inv.ef.rds")

# add a little variation to where the points in the violin plot are,
# makes it more visually appealing
set.seed(42)
map$dummy <- 1
map$dummy[map$Disease.state != "Patient" ] <- 2
map$dummy  <- map$dummy + sample(-100 : 100, nrow(map), replace = T)/2000

map.pat <- map[map$Disease.state == "Patient", ]
map.rel <- map[map$Disease.state == "Relative", ]



vioplot(map.pat$invsimp,
        map.rel$invsimp,
        drawRect = F, col=unique(map$col),
        names = c("Patients", "Relatives"),
        ylab = "invsimp")


segments(.7, mean(map.pat$invsimp),
         1.3, mean(map.pat$invsimp), lwd = 2)

segments(1.7, mean(map.rel$invsimp),
         2.3, mean(map.rel$invsimp), lwd = 2)

families <- unique(map$Family.ID)

for(i in 1: length(families)){

  map.pick <- map.pat[map.pat$Family.ID ==  families[i] , ]
  map.pick2 <- map.rel[map.pat$Family.ID ==  families[i] , ]

  segments(y0 = map.pick$invsimp, x0 = map.pick$dummy,
           y1 = map.pick2$invsimp, x1 = map.pick2$dummy,
           col = 8)

  points(map.pick$dummy, map.pick$invsimp, pch = 21, bg = 8)
  points(map.pick2$dummy, map.pick2$invsimp, pch = 21, bg = 8)


}
title(paste("Paired permutational t-test\nmean p value =",
            round(mean(t.inv.p, na.rm = T), 5), ", mean effect size = " ,
            round(mean(t.inv.ef, na.rm = T), 5),
            "\nmean 95% CI =", round(mean(t.inv.lo), 3),
            round(mean(t.inv.hi), 3), "
            df = ", mean(t.inv.df)))

#############################################################################
#               PANEL B distance to group centroid
#############################################################################


t.beta.stat <- t.beta.p <- t.beta.df <- t.beta.ef <- NULL

pats <- unique(map$DonorID[map$Disease.state == "Patient"])

#compare distance to group centroid
for(i in  1 : num.perm){

  print(i)

  map.i <- paired.map(map, pats)

  t.beta.res <- t.test( map.i$dist.cent[map.i$Disease.state == "Patient"],
                        map.i$dist.cent[map.i$Disease.state == "Relative"],
                        paired = T)

  t.beta.p[i] <-  t.beta.res$p.value
  t.beta.ef[i] <- t.beta.res$estimate
  t.beta.df[i] <- t.beta.res$parameter
  t.beta.stat[i] <- t.beta.res$statistic
}

dir.create("Microbiome_Permution_Analysis")

saveRDS(t.beta.p, "Microbiome_Permution_Analysis/t.beta.p.rds")
saveRDS(t.beta.df, "Microbiome_Permution_Analysis/t.beta.df.rds")
saveRDS(t.beta.ef, "Microbiome_Permution_Analysis/t.beta.ef.rds")
saveRDS(t.beta.stat, "Microbiome_Permution_Analysis/t.beta.stat.rds")


t.beta.p <- readRDS("Microbiome_Permution_Analysis/t.beta.p.rds")
t.beta.df <- readRDS("Microbiome_Permution_Analysis/t.beta.df.rds")
t.beta.ef <- readRDS("Microbiome_Permution_Analysis/t.beta.ef.rds")

set.seed(7)
map$dummy <- 1
map$dummy[map$Disease.state != "Patient" ] <- 2
map$dummy  <- map$dummy + sample(-100 : 100, nrow(map), replace = T)/2000

map.pat <- map[map$Disease.state == "Patient", ]
map.rel <- map[map$Disease.state == "Relative", ]


vioplot(map$dist.cent ~ map$Disease.state,
        drawRect = F, col=unique(map$col),
        names = c("Patients", "Relatives"),
        ylab = "BC distance to group centroid")


segments(.7, mean(map.pat$dist.cent),
         1.3, mean(map.pat$dist.cent), lwd = 2)

segments(1.7, mean(map.rel$dist.cent),
         2.3, mean(map.rel$dist.cent), lwd = 2)

families <- unique(map$Family.ID)

for(i in 1: length(families)){

  map.pick <- map.pat[map.pat$Family.ID ==  families[i] , ]
  map.pick2 <- map.rel[map.pat$Family.ID ==  families[i] , ]

  segments(y0 = map.pick$dist.cent, x0 = map.pick$dummy,
           y1 = map.pick2$dist.cent, x1 = map.pick2$dummy,
           col = 8)

  points(map.pick$dummy, map.pick$dist.cent, pch = 21, bg = 8)
  points(map.pick2$dummy, map.pick2$dist.cent, pch = 21, bg = 8)


}
title(paste("Paired T test mean p value =",
            round(mean(t.beta.p, na.rm = T), 5), "\ndf = ",
            paste(unique(t.beta.df), collapse = " "),
            "\nmean of the differences",
            round(mean(t.beta.ef, na.rm = T), 5)))

########## compile source data file ###########33
sup_data_1a <- cbind(map.pick$sample.id,
                     map.pick$DonorID,
                     map.pick$Family.ID,
                     map.pick$Disease.state,
                     map.pick$invsimp)

colnames(sup_data_1a) <- c("sample.id",
                           "donor.id",
                           "family.id",
                           "Disease State",
                           "Inverse Simpson Index")

sup_data_1a <- cbind("", "", sup_data_1a)
sup_data_1a[1,1] <- "Supplemental figure 1b"

write.csv(sup_data_1a,"sup_data_1a.csv", row.names = F)


sup_data_1b <- cbind(map.pick$sample.id,
                        map.pick$DonorID,
                        map.pick$Family.ID,
                        map.pick$Disease.state,
                        map.pick$dist.cent)

colnames(sup_data_1b) <- c("sample.id",
                              "donor.id",
                              "family.id",
                              "Disease State",
                              "Bray-Curtis distance to group centroid")

sup_data_1b <- cbind("", "", sup_data_1b)
sup_data_1b[1,1] <- "Supplemental figure 1b"

write.csv(sup_data_1b,"sup_data_1b.csv", row.names = F)





