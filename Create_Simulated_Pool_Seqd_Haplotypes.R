#### This file contains an R script described on the Step 3 of the ReadMe file
### This R script runs a single simulation of a "Pool-Seq" event.  Simulated cells that are sampled at the end of the experiment are input at Step 1 (of this script) and pooled
### At step 2 (of this script), cells are assinged a random haplotype from a pool of 60,000 haplotypes that was generated previously ("Create_All_Haplotype_Using_Site_Frequency_Spectrum.R")

args <- commandArgs(trailingOnly = TRUE)
## Derives two arguments from the slurm script that is queued, simulation ID and the R 'seed' 

seed_sim=args[1]
sim_id=args[2]
set.seed(seed_sim)

##### Step 1
##### Pooling cells to simulate Pool-Seq
##### cell_number_??? contains the abundance of each cell (from the original population) that is present at the end of the experiment
cell1 <- gsub("cellno", sample(1:200)[1], "$ADDRESS/cell_number_cellno")
cell2 <- gsub("cellno", sample(1:200)[1], "$ADDRESS/cell_number_cellno")
cell3 <- gsub("cellno", sample(1:200)[1], "$ADDRESS/cell_number_cellno")

cellno1 <- read.table(cell1)
cellno2 <- read.table(cell2)
cellno3 <- read.table(cell3)

pool <- cellno1 + cellno2 + cellno3


rm(cellno1); rm(cellno2); rm(cellno3);

## This file contains an R script described on the Step 3 of the ReadMe file
## "Master.matrix.Ancestor" contains 60,000 haplotypes that are consistent with the site frequency spectrum observed in the ancestral population


cAnc <- read.csv("$Master.matrix.Ancestor")

## Pre-allocating a matrix to contain 60,000 haplotypes {2*sum(pool)} across 42004 loci {length(cAnc[sample(1:nrow(cAnc))[1], ])))}
master_matrix <- data.frame(matrix(nrow=(2*sum(pool)), ncol=length(cAnc[sample(1:nrow(cAnc))[1], ])))


count=0;
Start <- Sys.time();
for (i in 1:nrow(pool)){
#TESTfor (i in 1:200){
realized.haplotype_1 <- cAnc[sample(1:nrow(cAnc))[1], ]
realized.haplotype_2 <- cAnc[sample(1:nrow(cAnc))[1], ]


        if (pool$x[i]==0) {next}
        else for (j in 1:pool$x[i]){
        count=count+1;

        master_matrix[count,] <- realized.haplotype_1
        count=count+1; ## Accomodates the diploid haplotype2
        master_matrix[count,] <- realized.haplotype_2


}

print(sprintf("Total time so far processing %s cells", i))
print(Sys.time()-Start)

}

write.csv(master_matrix, gsub("sim_id", sim_id, "$Pooled Sequences . csv"))

## This csv file contains 60,000 haplotypes, from 30,000 cells after pooling 3 evolved populations




