The Simulations were run using four scripts, in this order

1.  "Generate_Neutral_Distribution_Of_Cells_After_90_days.txt" is an R script that models the neutral evolution of 10,000 cells,
through doubling events and bottlenecks, and outputs a frequency distribution of all 10,000 cells after the experiment ends in 90 days.
Roughly 1800 cells survive at the end of each simulation, and the rest are lost to genetic drift (random sampling in a finite population).
This script was run 1,000 times to generate a thousand instances of neutral evolution.

2.  "Create_All_Haplotype_Using_Site_Frequency_Spectrum.txt" is an R script that uses the empirically observed site frequency spectrum
(SFS; Allele frequency w.r.t loci) and generates a population of 60,000 haplotypes that are consistent with the observed SFS.
This is an approximation of the pool of ancestral cells before the experiment began.  Because the linkage structure of the original population 
is unknown, the haplotypes are assigned genotypes at each loci independently (leading to near-complete linkage equilibrium)

3.  "Create_Simulated_Pool_Seqd_Haplotypes.txt" is an R script that takes the input from the previously run scripts to simulate "Pool-Seq".  
It chooses 3 instances of neutral evolution from Step 1, pools them together, and assigns each cell a genotype using the haplotypes generated 
by the second script.  The genotypes could have been assigned at any stage of the simulation as all cells are considerered phenotypically identical

4.  "Add_Sequencing_Variation.txt"  is an R script that analyzes the Pool-Seq'ed samples after adding sequencing variation to them.  
This generates the final simulated allele frequencies, p*, by locus.