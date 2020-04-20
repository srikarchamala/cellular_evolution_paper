#### This file contains an R script described on the Step 4 of the ReadMe file
#### Once the pooled haplotypes are created, we can estimate final allele frequencies by simply averaging them over (by each loci)
#### If we sample, out of the 60,000 haplotypes, only ~250x (equivalent to the coverage at the loci in the evolved lines), 
#### then we account for sequencing variation at the level of the evolved lines (generating p')
#### The Binomial sampling step, generates p* from p' by using the read coverage at the loci in the ancestor (described in the Methods section)


n <- length(list_of_pooled_haplotypes_sets)
## n=60000



allele_freq <- data.frame(matrix(nrow=n, ncol=42004))
allele_freq_sample_evolved <- data.frame(matrix(nrow=n, ncol=42004)) ## 42004 is the number of quality filtered variants
allele_freq_sample_evolved_Anc <- data.frame(matrix(nrow=n, ncol=42004))

cov_data <- read.csv("Coverage data of all lines", sep="\t", header=TRUE)
cov_evolved <- as.vector((cov_data)$Ave_Soft_Stiff)
cov_ancestor <- as.vector((cov_data)$Anc)
#Start <- Sys.time()
for (i in 1:n){
        print(list_of_pooled_haplotypes_sets[i])
        pooled_haplotype_temp <- read.csv(list_of_pooled_haplotypes_sets[i])
        pooled_haplotypes  <- pooled_haplotype_temp[,-1] ; rm(pooled_haplotype_temp)

        print(nrow(pooled_haplotypes))
        pooled_haplotypes_temp <- pooled_haplotypes[complete.cases(pooled_haplotypes[c(1:nrow(pooled_haplotypes)),]),]
        pooled_haplotypes <- pooled_haplotypes_temp; rm(pooled_haplotypes_temp)
        print(nrow(pooled_haplotypes))

        #theta_ar[i] <- theta(pooled_haplotypes)
        #pie_ar[i] <- pie2(pooled_haplotypes)

        ## Af change with subsampling of reads

        for (col in 1:42004){
		allele_freq_sample_evolved[i,col] <- mean(sample(pooled_haplotypes[,col], cov_evolved[col], replace=TRUE), na.rm=TRUE)
        samp.1 <- mean(sample(pooled_haplotypes[,col], cov_evolved[col], replace=TRUE), na.rm=TRUE)
        samp.2 <- mean(sample(rbinom(60000,1,samp.1),cov_ancestor[col], replace=TRUE))
        allele_freq_sample_evolved_Anc[i,col] <- samp.2

        }

        allele_freq[i,] <- colMeans(pooled_haplotypes, na.rm = FALSE, dims = 1)

        print(i)


}
write.csv(allele_freq, gsub("total_simulations", n, "$Address/Allele_frequency_change_table_NO_SEQ_VAR_neutral_pooled_evolution.total_simulations.csv"))
write.csv(allele_freq_sample_evolved, gsub("total_simulations", n, "$Address/Allele_frequency_change_table_EVOLVED_neutral_pooled_evolution.total_simulations.csv"))
write.csv(allele_freq_sample_evolved_Anc, gsub("total_simulations", n, "$Address/Allele_frequency_change_table_ANC_EVOLVED_neutral_pooled_evolution.total_simulations.csv"))

## We generated three files, each containing the final allele frequency
## First with no sequencing variation, then with sequencing variation considered only in the evolved lines, and thirdly sequencing variation considered both in the evolved lines and the ancestor.