#### Step 2 from the ReadMe file
#### This script creates a population of haplotypes given a tab delimited file with allele frequencies and uses a vcf file to filter sites by quality
#### This is an R script that is run on a SLURM scheduling system using the following command -
#### R CMD BATCH --vanilla "This following script"



#### Template delimited file
#chr     pos     ref     alt     SOFT2   SOFT1   STIFF1  ANCESTRAL
#1       3206823 G       T       0.5255131964809384      0.5277401894451962      0.6663101604278074      0.5930630190522717
#1       3206889 T       G       0.011547344110854504    0.007069408740359897    0.014361702127659574    0.008687258687258687
#1       3411768 A       C       0.008964955175224124    0.00872093023255814     0.014221556886227544    0.009272467902995721
#1       3411930 C       T       0.40822590546347454     0.4455713319810683      0.5595026642984015      0.5185766614338043
#1       3660709 T       G       0.010572687224669603    0.0009624639076034649   0.006704980842911878    0.0029239766081871343
#1       3661014 T       G       0.010507880910683012    0.013295346628679962    0.00427715996578272     0.018138801261829655


library(dplyr)


data <- read.delim("$Tab_delimited_file")
## Data curation
data$ANCESTRAL <- as.numeric(as.character(data$ANCESTRAL)); data$ANCESTRAL[is.nan(data$ANCESTRAL)] <- NA  ;
data$SOFT2 <- as.numeric(as.character(data$SOFT2)); data$SOFT2[is.nan(data$SOFT2)] <- NA  ;
data$SOFT1 <- as.numeric(as.character(data$SOFT1)); data$SOFT1[is.nan(data$SOFT1)] <- NA  ;
data$STIFF1 <- as.numeric(as.character(data$STIFF1)); data$STIFF1[is.nan(data$STIFF1)] <- NA  ;


## This step filters out any loci that consistently had >99% or <1% allele frequency
data <- subset(data, data$SOFT1 < 0.99 & data$STIFF1 < 0.99 & data$SOFT2 < 0.99 & data$ANCESTRAL < 0.99)
data <- subset(data, data$SOFT1 > 0.01 & data$STIFF1 > 0.01 & data$SOFT2 > 0.01 & data$ANCESTRAL > 0.01)

##vcf filters
vcf_data <- read.delim("$vcf.R FILE") ## The .R file has all the "#" fields removed so that a vcf file can be read like a table
vcf_data1 <- subset(vcf_data, QUAL > 30)

colnames(data) <- c("X.CHROM", "POS", "REF", "ALT", "soft2", "soft1", "stiff1", "ancestral") ## making the soft1/2 colnames different from the vcf before inner join so that they don't interfere
data <- semi_join(data, vcf_data1)
colnames(data) <- c("X.CHROM", "POS", "REF", "ALT", "SOFT2", "SOFT1", "STIFF1", "ANCESTRAL") ## Getting column names back to normal

### Remove Ancestor with NA's
nrow(data)
data.Anc <- data[complete.cases(data[,c(5,6,7,8)]),]

## In this section, we will create a master matrix that contains the genotype of 30,000 cells or 60,000 haplotypes

nSites=length(data.Anc$ANCESTRAL)
sfs=data.Anc$ANCESTRAL
nCells=30000
realized_genotype_haploid <- array()
master_matrix <- rep(0,nSites)


for (cell_ID in 1:(2*nCells)){

  for (i in 1:nSites){
    realized_genotype_haploid[i] <- ifelse(runif(1) <= sfs[i], 1, 0)
  }
print(cell_ID);

  assign(gsub("Haploid_ID",cell_ID, "realized_genotype_haploid_Haploid_ID"), realized_genotype_haploid)
  master_matrix <- rbind(master_matrix,get(gsub("Haploid_ID",cell_ID, "realized_genotype_haploid_Haploid_ID")))
}

master_matrix <- master_matrix[-1,] ## Removing the blank row
write.csv(master_matrix, gsub("sim", j, "$Master_matrix_Saved_as_csv"))
}

