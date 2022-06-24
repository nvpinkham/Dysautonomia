# Dysautonomia
R code used in "Gut microbiome dysbiosis drives metabolic dysfunction in Familial dysautonomia"

Each sample's meta data is found in the "data" folder alongside metabolite and OTU tables
Raw sequencing data is avilible at under bioproject PRJNA785599 https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA785599

R folder contains the custom R
Man folder contains function documentation 

The scripts for creating each figure are in the associated folders
 - these scripts are supposed to be run in their parent directory (same directory as this document)

### System requirements
This analysis was done with macOS version 12.4
R version 4.1.0

MetaboAnalyst was not easily compiled on our system so we used individual R scripts from the package. Available here https://github.com/xia-lab/MetaboAnalystR 

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
 
