
venn2otu <- function(otu, map,
                     group1 = "mouse gut metagenome",
                     group2 = "human gut metagenome",
                     min.counts = 0,
                     min.sites = 5,
                     out.file="ExampleVenn.EPS",
                     export = F){

  if(!"split" %in% colnames(map)){
    message('need col named "split"')
  }

  if(!all(row.names(map) == row.names(otu))){
    return('CHECK that samples match')
  }

  otu.g1 <- otu[map$split == group1, ]
  otu.g2 <- otu[map$split == group2, ]

  aa <- min(c(nrow(otu.g1),
              nrow(otu.g2)))

  #set.seed(42)
  otu.g1 <- otu.g1[sample(1:nrow(otu.g1), aa),]
  otu.g2 <- otu.g2[sample(1:nrow(otu.g2), aa),]

  g1 <- colSums(otu.g1)
  g2 <- colSums(otu.g2)

  g1.pa <- colSums(otu.g1 > min.counts)
  g2.pa <- colSums(otu.g2 > min.counts)

  unique2g1 <- g1.pa > min.sites & g2.pa < min.sites
  unique2g2 <- g2.pa > min.sites & g1.pa < min.sites

  shared <- g1.pa >= min.sites & g2.pa >= min.sites

  sum(otu[ , shared])
  sum(otu[ , unique2g1])
  sum(otu[ , unique2g2])

  if(export){

    setEPS()
    postscript(file = out.file, fonts = "serif")
    #grid::grid.newpage()
    try(venn.plot <- VennDiagram::draw.pairwise.venn(area1=length(unique2g1[unique2g1]),
                                                     area2=length(unique2g2[unique2g2]),
                                                     cross.area=length(shared[ shared]),
                                                     category=c(group1, group2),
                                                     fill=c("Red","Yellow")))

    try(grid::grid.draw(venn.plot))

    dev.off()

  }else{

    grid::grid.newpage()

    venn.plot <- VennDiagram::draw.pairwise.venn(area1=length(unique2g1[unique2g1]),
                                                 area2=length(unique2g2[unique2g2]),
                                                 cross.area=length(shared[ shared]),
                                                 category=c(group1, group2),
                                                 fill=c("Red","Yellow"))

    grid::grid.draw(venn.plot)
  }

  g1.proportion <- sum(otu.g1[ , shared]) /sum(otu.g1)
  g2.proportion <- sum(otu.g2[ , shared]) /sum(otu.g2)

  res <- c(g1.proportion, g2.proportion)
  names(res) <- c(group1, group2)
  print(res)


  # res <- list(tax[unique2g1,], tax[shared,],tax[unique2g2,])

  # names(res) <- c(group1,"shared", group2)

  return(res)
}




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




