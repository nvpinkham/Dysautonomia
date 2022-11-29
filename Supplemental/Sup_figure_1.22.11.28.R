

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


t.inv <- as.data.frame(matrix(nrow = num.perm, ncol = 6))
colnames(t.inv) <- c("P_val", "df", "effect_size", "t.stat", "conf_int_lo",  "conf_int_hi")

pats <- unique(map$DonorID[map$Disease.state == "Patient"])

for(i in  1 : num.perm){
  # paired permutation analysis
  # Each iteration of loop is a permutation

  print(i)

  map.i <- paired.map(map, pats) # randomly sampled each time

  res <- t.test( map.i$invsimp[map.i$Disease.state == "Patient"],
                       map.i$invsimp[map.i$Disease.state == "Relative"],
                       paired = T)

  t.inv$P_val[i] <- res$p.value
  t.inv$df[i] <- res$parameter
  t.inv$effect_size[i] <- res$estimate
  t.inv$t.stat[i] <- res$statistic
  t.inv$conf_int_lo[i] <- res$conf.int[1]
  t.inv$conf_int_hi[i] <- res$conf.int[2]
}

# save results from permutation analysis
dir.create("Microbiome_Permution_Analysis")

saveRDS(t.inv, "Microbiome_Permution_Analysis/t.inv.rds")

t.inv <- readRDS("Microbiome_Permution_Analysis/t.inv.rds")

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

title(paste("Paired T test mean p value =",  round(mean(t.inv$P_val, na.rm = T), 5),
            "\nmean effect size =", round(mean(t.inv$effect_size, na.rm = T), 5),
            "\ndf =", round(mean(t.inv$df, na.rm = T), 5),
            "\nmean 95% CI =", round(mean(t.inv$conf_int_lo, na.rm = T), 2),"to",
            round(mean(t.inv$conf_int_hi, na.rm = T), 2)))

#############################################################################
#               PANEL B distance to group centroid
#############################################################################
beta.t <- as.data.frame(matrix(nrow = num.perm, ncol = 6))
colnames(beta.t) <- c("P_val", "df", "effect_size", "t.stat", "conf_int_lo",  "conf_int_hi")



pats <- unique(map$DonorID[map$Disease.state == "Patient"])

#compare distance to group centroid
for(i in  1 : num.perm){

  print(i)

  map.i <- paired.map(map, pats)

  res <- t.test( map.i$dist.cent[map.i$Disease.state == "Patient"],
                 map.i$dist.cent[map.i$Disease.state == "Relative"],
                 paired = T)

  beta.t$P_val[i] <- res$p.value
  beta.t$df[i] <- res$parameter
  beta.t$effect_size[i] <- res$estimate
  beta.t$t.stat[i] <- res$statistic
  beta.t$conf_int_lo[i] <- res$conf.int[1]
  beta.t$conf_int_hi[i] <- res$conf.int[2]
}

dir.create("Microbiome_Permution_Analysis")

saveRDS(beta.t, "Microbiome_Permution_Analysis/t.beta.rds")

beta.t <- readRDS("Microbiome_Permution_Analysis/t.beta.rds")

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

title(paste("Paired T test mean p value =",  round(mean(beta.t$P_val, na.rm = T), 5),
            "\nmean effect size =", round(mean(beta.t$effect_size, na.rm = T), 5),
            "\ndf =", round(mean(beta.t$df, na.rm = T), 5),
            "\nmean 95% CI =", round(mean(beta.t$conf_int_lo, na.rm = T), 2),"to",
            round(mean(beta.t$conf_int_hi, na.rm = T), 2)))

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


t.inv.means <- apply(t.inv, 2, mean)
names(t.inv.means) <- paste("mean over 999 permutations", names(t.inv.means))
t.inv.means <- capture.output(t.inv.means)
writeLines(t.inv.means, "Sup_Fig_1a_tests.txt")

beta.t.means <- apply(beta.t, 2, mean)
names(beta.t.means) <- paste("mean over 999 permutations", names(beta.t.means))
beta.t.means <- capture.output(beta.t.means)
writeLines(beta.t.means, "Sup_Fig_1b_tests.txt")




