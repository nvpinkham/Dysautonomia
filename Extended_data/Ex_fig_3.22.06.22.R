

setwd("/Users/nickpinkham/Desktop/Nature2022_submission/Extended_Data/")

source("/Users/nickpinkham/Desktop/Nature2022_submission/Figure_1/human_functions.22.05.04.R")

## PAIRED COMPARISONS FROM HERE :

map$Collection.date <- as.Date(map$Collection.date, format = "%m/%d/%y")
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


res.tax <- tax.shared(map, otu, tax)
dim(res.tax)

fams <- unique(tax$family)

res.array <- array(dim = c(length(fams), 120, 999))

rownames(res.array) <- rownames(res.tax)
colnames(res.array) <- colnames(res.tax)
res.array.discrete <- res.array

for(i in 1:999){
  print(i)
  res.array[ , , i] <- tax.shared(map, otu, tax)
  res.array.discrete[ , , i] <-  tax.shared(map, otu, tax, binary = T)
}

dir.create("permution_analysis")
saveRDS(res.array, "Family_array.rds")
saveRDS(res.array.discrete, "Family_array_discrete.rds")

res.array <- readRDS("Family_array.rds")
res.array.discrete <- readRDS("Family_array_discrete.rds")

res.m <- rowMeans(res.array, dims = 2, na.rm = T)
res <- aggregate(t(res.m),
                 list(colnames(res.m)), 
                 mean)

rownames(res) <- res$Group.1
res$Group.1 <- NULL

res.log <- log10(as.matrix(res) + 1)#+ .000001)
res.log <- res.log[ , order(colSums(res.log)) ]


res.array <- readRDS("Family_array.rds")
res.array.pa <- readRDS("Family_array_discrete.rds")

res.tax <- rowMeans(res.array, dims = 2, na.rm = T)
res.tax.pa <- rowMeans(res.array.pa, dims = 2, na.rm = T)


#res.tax <- tax.shared(map, otu, tax)
res.tax <- aggregate(t(res.tax),
                     list(colnames(res.tax)), 
                     mean)

row.names(res.tax) <- res.tax$Group.1
res.tax$Group.1 <- NULL

res.log <- log10(as.matrix(res.tax) + 1)#+ .000001)
res.log <- res.log[  , order(colSums(res.log)) ]

res.tax.pa <- aggregate(t(res.tax.pa),
                        list(colnames(res.tax.pa)), 
                        mean)

row.names(res.tax.pa) <- res.tax.pa$Group.1
res.tax.pa$Group.1 <- NULL
res.tax.pa <- res.tax.pa[-2,]
res.tax.pa <- res.tax.pa[ , match(colnames(res.log), 
                                  colnames(res.tax.pa))]

row.names(res.tax.pa) <- c("  P", "  S", "  R")



par(mar = c(3, 15, 4, 4))

par(xpd = TRUE) #Draw outside plot area

barplot(res.log, 
         xlim = c(-3.3, 16),
         cex.axis=.75, 
         cex.names=.75,
         xaxt='n',
         horiz=T ,
         main = "FAMILY",
         xlab = "median population size (log 10)",
         las = 1,
         col =c( "#E69F00","grey65", "gray77", "#56B4E9"))

# text(-4.5, 1:42 * 1.2 - .5, round(rev(ps), 5))

row.names(res.tax.pa) <- c(" Patient", " Shared", " Relative")

  col.fun <- heat(matrix = res.tax.pa, 
                  colors =  c("purple", "yellow", "orange", "orangered"),
                  x.adj =  -3.25)
  
  labs <- round(seq(min(res.tax.pa), max(res.tax.pa),l=5))
  text(x=6.5, y = seq(5,10,l=5), labels = labs, cex = .5)
  text(5, 10.5, "# OTUs", cex = 0.5)
  
  legend_image <- as.raster(matrix(col.fun(100), ncol=1))
  
  rasterImage(rev(legend_image),
              xleft = 4, 
              xright = 6,
              ytop =  10, 
              ybottom = 5)
  
  
  segments(x0 = 12, x1 = 13, y0 = 14, y1 = 14)
  text(12.5, 15, "10", cex = 0.5)
  
  segments(x0 = 11, x1 = 13, y0 = 12, y1 = 12)
  text(12, 13, "100", cex = 0.5)
  
  segments(x0 = 10, x1 = 13, y0 = 10, y1 = 10)
  text(11.5, 11, "1000", cex = 0.5)
  
  op <- par(cex = 0.5)
  
  legend("bottomright", 
         fill = c( "#E69F00", "grey65", "gray78", "#56B4E9"),
         bty = "n",
         legend = rev(c("Unique in Relative", 
                        "Shared within pair; in Relative", 
                        "Shared within pair; in Patient",
                        "Unique in Patient")))
