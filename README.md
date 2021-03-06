# Dysautonomia
The content of this GitHub entry was used in the manuscript, "Gut microbiome dysbiosis drives metabolic dysfunction in Familial dysautonomia"

1) The "R" folder contains the custom scripts for analyzing human and mouse experiments

2) The "man" folder contains function documentation 

3) The "data" folder contains metadata along with normalized metabolite abundances and 16S rRNA sequencing OTU tables

4) Folders labeled "Figure_" contain R scripts for generating each figure shown in the main manuscript. Run script in the parent directory.
 - these scripts are supposed to be run in their parent directory (same directory as this document)

Raw sequencing data is available under bioproject PRJNA785599 https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA785599 
Raw metabolomic data is available MetaboLights accession TBD for now the raw NMR data is available on OneDrive:  it is password protected 

### System requirements
This analysis was done with macOS version 12.4
R version 4.1.0
MetaboAnalyst is not easily compiled on current macOS so individual R scripts from this package were used for certain processes. These are available at https://github.com/xia-lab/MetaboAnalystR 

if you are having trouble with the MetaboAnalyst scripts we have noticed that the R package cairo is often the cause. It is not needed for any of the functions from MetaboAnalyst we utilized so the lines in MetaboAnalyst's "generl_data_utils" that call cairo can be deleted to get around this problem. 

### Dependancies
  randomForest 
  rfUtilities
  vegan 
  vegan3d
  labdsv
  cluster
  viopot
  RColorBrewer
  dendextend
  lme4
  mixOmics
  dplyr

### Time estimates for R scripts associated with each figure:

 figure 1:\
   a. < 10 sec.\
  b. 10 - 15 minutes for permutation analysis * change number "num.perm" to run few permutations to speed up process\
  c. 10 - 15 minutes for permutation analysis\
  d. < 10 sec.\
  e. < 10 sec.
  
 figure 2:\
  a. 3 - 5 minutes for Random Forest annalysis\
  b. 20 - 30 minutes for permutation analysis\
  c. < 10 sec.\
  d. < 10 sec.
  
 figure 3:\
  a. < 10 sec.\
  b. < 10 sec.\
  c. < 10 sec.\
  d. < 10 sec.
  
 figure 4:\
  a. 5 minutes for multiple PERMANOVAs \
  b. 10 minutes for multiple PERMANOVAs \
  c. < 10 sec.\
  d. < 10 sec.
 
### mothur pipeline 
The R scripts used to make the human and mice OTU tables are included in the "Run_mothur" folder



