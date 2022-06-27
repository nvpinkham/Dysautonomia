

getwd()

dir.create("processing")

source("/Users/nickpinkham/Desktop/Nature2021_final/Figure_2/Fig_2_functions.11.10.R")

map$illumina.Sample.Alias

reads.path <- list.files("/Users/nickpinkham/Desktop/Nature2021_final/data/human/Illumina_data/Raw_Human_Reads", full.names = T)
reads <- list.files("/Users/nickpinkham/Desktop/Nature2021_final/data/human/Illumina_data/Raw_Human_Reads")


ids <- gsub("_S.*_L001_R.*_001.fastq.gz", "", reads)

ids <- gsub("[-]|[.]", "_", ids)

map$illumina.id <- gsub("-|[.]", "_", map$illumina.Sample.Alias)

map$illumina.Sample.Alias[i] 

for(i in 1 :  nrow(map)){
  
  map$r1[i] <- reads.path[ids == map$illumina.id[i]][1]
  map$r2[i] <- reads.path[ids == map$illumina.id[i]][2]
}


dir.create("/Users/nickpinkham/Desktop/Nature2021_final/data/human/Illumina_data/Final_Reads")
final_folder <- "/Users/nickpinkham/Desktop/Nature2021_final/data/human/Illumina_data/Final_Reads/"

for(i in 1 : nrow(map)){
  
  print(i)
  
  map$r1.final <- paste0(final_folder, map$sample.id[i], "_S", i, "_L001_R1.fastq.gz")
  map$r2.final <- paste0(final_folder, map$sample.id[i], "_S", i, "_L001_R2.fastq.gz")
  
  file.copy(from = map$r1[i], to = map$r1.final[i])
  file.copy(from = map$r2[i], to = map$r2.final[i])
}

################################################################
system('./mothur "#set.seed(seed = 42);
                   make.file(inputdir = human/Illumina_data/Final_Reads, type=gz, prefix=FD_human21);
                   make.contigs(file=current, processors = 12);
                   rename.seqs(fasta=current, group=current);"')

  file.copy(from = "human/Illumina_data/Final_Reads/FD_human21.trim.contigs.renamed.fasta", 
            to = "processing/FD_human21.trim.contigs.renamed.fasta")# move to processing directory
  
  file.copy(from = "human/Illumina_data/Final_Reads/FD_human21.contigs.renamed.groups", 
            to = "processing/FD_human21.contigs.renamed.groups")      
  
system('./mothur "#set.dir(input=processing);
                   screen.seqs(fasta = FD_human21.trim.contigs.renamed.fasta, maxambig=0, maxlength=253, maxhomop=8);
                   unique.seqs();
                   count.seqs(name=current, group = FD_human21.contigs.renamed.groups);
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
  
file.copy(from = "processing/FD_human21.trim.contigs.renamed.good.unique.good.filter.precluster.unique.pick.pick.agc.shared", 
          to = "human/FD_human21.shared")   

file.copy(from = "processing/FD_human21.trim.contigs.renamed.good.unique.good.filter.precluster.unique.pick.pick.agc.0.03.cons.taxonomy", 
          to = "human/FD_human21.taxonomy")  

tax <- read.table("human/FD_human21.taxonomy", row.names = 1, header = T)

otu <- read.table("human/FD_human21.shared", 
                  row.names = 2, header = T)
otu$label <- NULL
otu$numOtus <- NULL

tax <- tax[   colSums(otu) > 100 , ]
otu <- otu[ , colSums(otu) > 100   ]

min(rowSums(otu))
otu <- vegan::rrarefy(otu, 9848) 

write.csv(otu, "human/FD_OTU_human21.csv")
  
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

write.csv(tax.df, "human/FD_OTU_human21_taxonomy.csv")
################################################################################


