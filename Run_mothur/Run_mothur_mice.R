

getwd()
dir.create("processing")

#source("/Users/nickpinkham/Desktop/Nature2021_final/Figure_4/fig4_functions.11.23.21.R")
#map <- map[!is.na(map$X16SFile.ID) , ]


map.mice.3 <- read.csv("mice/ASD_map.csv",
                       header=TRUE, row.names=1)

map.mice.3$sample.id <- paste0(map.mice.3$mouse, 
                                 ".",
                                 substr(map.mice.3$genotype, 1, 1),
                                 ".",
                                 map.mice.3$time)

map.mice.4 <- read.csv("mice/CoHoused_mice.8.16.csv",
                       header=TRUE)
map.mice.4$sample.id <- paste0("X", map.mice.4$ID.Gen.Tmt.Age)

map.mice.4 <- map.mice.4[!is.na(map.mice.4$X16SFile.ID) & 
                           map.mice.4$X16SFile.ID != "" , ]

reads.path <- list.files("/Users/nickpinkham/Desktop/Nature2021_final/data/mice/Illumina_data/Raw_Mouse_Reads", 
                         full.names = T, pattern = "fastq")
reads <- list.files("/Users/nickpinkham/Desktop/Nature2021_final/data/mice/Illumina_data/Raw_Mouse_Reads",
                     pattern = "fastq")

ids <- gsub("_S.*_L001_R.*_001.fastq.gz", "", reads)
ids <- gsub("[-]|[.]", "_", ids)

map.mice.4$illumina.id <- gsub("-|[.]", "_", map.mice.4$X16SFile.ID)

for(i in 1 :  nrow(map.mice.4)){
  
  map.mice.4$r1[i] <- reads.path[ids == map.mice.4$illumina.id[i]][1]
  map.mice.4$r2[i] <- reads.path[ids == map.mice.4$illumina.id[i]][2]
}

r1 <- c(paste0("/Users/nickpinkham/Desktop/Nature2021_final/data/mice/Illumina_data/Raw_Mouse_Reads/", map.mice.3$r1), 
        map.mice.4$r1)
r2 <- c(paste0("/Users/nickpinkham/Desktop/Nature2021_final/data/mice/Illumina_data/Raw_Mouse_Reads/", map.mice.3$r2), 
        map.mice.4$r2)

sample.id <- c(map.mice.3$sample.id, map.mice.4$sample.id)

map <- as.data.frame(cbind(sample.id, r1, r2))

dir.create("/Users/nickpinkham/Desktop/Nature2021_final/data/mice/Illumina_data/Final_Reads")
final_folder <- "/Users/nickpinkham/Desktop/Nature2021_final/data/mice/Illumina_data/Final_Reads/"

gz <- substr(map$r1, nchar(map$r1) - 1, nchar(map$r1)) == "tq"
map$r1[gz] <- paste0(map$r1[gz], ".gz")
map$r2[gz] <- paste0(map$r2[gz], ".gz")

for(i in 1 : nrow(map)){
  
  print(i)
  
  map$r1.final[i] <- paste0(final_folder, map$sample.id[i], "_S", i, "_L001_R1.fastq.gz")
  map$r2.final[i] <- paste0(final_folder, map$sample.id[i], "_S", i, "_L001_R2.fastq.gz")

  file.copy(from = map$r1[i], to = map$r1.final[i])
  file.copy(from = map$r2[i], to = map$r2.final[i])

  map.mice.3$r1.final[map.mice.3$sample.id == map$sample.id[i]] <- map$r1.final[i]
  map.mice.3$r2.final[map.mice.3$sample.id == map$sample.id[i]] <- map$r2.final[i]
  
  map.mice.4$r1.final[map.mice.4$sample.id == map$sample.id[i]] <- map$r1.final[i]
  map.mice.4$r2.final[map.mice.4$sample.id == map$sample.id[i]] <- map$r2.final[i]
  
  
}

write.csv(map.mice.3, "mice/ASD_mice.12.07.csv")
write.csv(map.mice.4, "mice/CoHoused_mice.12.07.csv")


################################################################
system('./mothur "#set.seed(seed = 42);
                   make.file(inputdir = mice/Illumina_data/Final_Reads, type=gz, prefix=FD_mice21);
                   make.contigs(file=current, processors = 12);
                   rename.seqs(fasta=current, group=current);"')

file.copy(from = "mice/Illumina_data/Final_Reads/FD_mice21.trim.contigs.renamed.fasta", 
          to = "processing/FD_mice21.trim.contigs.renamed.fasta", overwrite = T)# move to processing directory

file.copy(from = "mice/Illumina_data/Final_Reads/FD_mice21.contigs.renamed.groups", 
          to = "processing/FD_mice21.contigs.renamed.groups", overwrite = T)      

system('./mothur "#set.dir(input=mice/Illumina_data/processing);
                   set.seed(seed = 42);
                   screen.seqs(fasta = FD_mice21.trim.contigs.renamed.fasta, maxambig=0, maxlength=253, maxhomop=8);
                   unique.seqs();
                   count.seqs(name=current, group = FD_mice21.contigs.renamed.groups);
                   align.seqs(fasta=current, reference=silva.v4.fasta);
                   summary.seqs();
                   screen.seqs(fasta=current, count=current, start=1968, end=11550, maxhomop=8);
                   filter.seqs(fasta=current, vertical=T, trump=.);
                   pre.cluster(fasta=current, count=current, diffs=3);
                   unique.seqs(fasta=current, count=current);
                   chimera.vsearch(fasta=current, count=current, dereplicate=t);
                   remove.seqs(fasta=current, accnos=current);
                   classify.seqs(fasta=current, count=current, reference=trainset16_022016.rdp.fasta, taxonomy=trainset16_022016.rdp.tax, cutoff=80);
                   remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);
                   dist.seqs(fasta=current, cutoff=0.03);
                   cluster(fasta=current, count=current, method=agc, cutoff=0.03);
                   make.shared(list=current, count=current, label=0.03);
                   classify.otu(list=current, count=current, taxonomy=current, label=0.03);"')
#   make.shared(list=current, count=current, label=1);
#  classify.otu(list=current, count=current, taxonomy=current, label=1);"')

file.copy(from = "processing/FD_mice21.trim.contigs.renamed.good.unique.good.filter.precluster.unique.pick.pick.agc.shared", 
          to = "mice/FD_mice21.shared", overwrite = T)   

file.copy(from = "processing/FD_mice21.trim.contigs.renamed.good.unique.good.filter.precluster.unique.pick.pick.agc.0.03.cons.taxonomy", 
          to = "mice/FD_mice21.taxonomy", overwrite = T)  

tax <- read.table("mice/FD_mice21.taxonomy", row.names = 1, header = T)

otu <- read.table("mice/FD_mice21.shared", 
                  row.names = 2, header = T)
otu$label <- NULL
otu$numOtus <- NULL

tax <- tax[   colSums(otu) > 100 , ]
otu <- otu[ , colSums(otu) > 100   ]

min(rowSums(otu))
set.seed(42)
otu <- vegan::rrarefy(otu, 5308) 

write.csv(otu, "mice/FD_OTU_mice21.csv")

tax.list <- strsplit(tax$Taxonomy, ";|[(]|[)]")

tax.df <- matrix(nrow = nrow(tax), ncol = 18)

for(i in 1:18){
  
  tax.df[ , i] <- sapply(tax.list , "[[", i)
  
}

tax.df <- as.data.frame(tax.df)
tax.df <- tax.df[ , tax.df[1,] != ""]

colnames(tax.df) <- c("Kingdom", "King.conf", 
                      "Phylum",  "Phy.conf",
                      "Class",   "Class.conf",
                      "Order",   "Order.conf", 
                      "family",  "fam.conf",
                      "genus",   "gen.conf")

write.csv(tax.df, "mice/FD_OTU_mice21_taxonomy.csv")
################################################################################
