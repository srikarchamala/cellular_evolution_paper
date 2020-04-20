## This file contains an R script described on the Step 1 of the ReadMe file
## The R script is run on the Slurm scheduler using the following command - module load R
R CMD BATCH --vanilla "R_script"


## A 'week' is 3-4 days worth of three doublings, which takes the population upto ~80k or so.  ## Historical artifact
## Initialize arrays
remaining_cells_count <- array()
remaining_cells_by_week <- array()
remaining_cells_by_week_by_simulations <- rep(0,26)
remaining_cell_lineages <- array()
remaining_cell_identity <- list()

for (simulations in 1:10000){
  print("Simulation without Ne")
  cell_number <- rep(1,10000)

for (weeks in 1:26){
  ## 26 bottlenecks in 13 weeks

  start.time <- Sys.time()
  cell.number.start.of.week <- cell_number
  # Modelling cell division for three rounds; when population reaches ~80,000 the bottleneck occurs

  for (cell_division in 1:3){cell_number = cell_number + rbinom(length(cell_number), cell_number, prob=1)}


print("Sum of cell numbers after 3 rounds of replication")
print(sum(cell_number))


######################################################################## Bottleneck ##########################

cell_fraction <- cell_number/sum(cell_number)
print("Sum of cell fraction, should be one")
print(sum(cell_fraction))


cell_number_def <- rep(0,10000)
print(head(c(which(!cell_number %in% c(0))),25));
for (bottleneck_population in 1:10000){
j=0;

## Using the 'runif' method, we add number of cells back in the new pool stochatically. To fill in the gap if 10k aren't reached, we randomly sample out 'remaining cells'
remaining_cells_shuffled <- sample(c(which(!cell_number %in% c(0))))
for (i in remaining_cells_shuffled){

r=runif(1,0,1); j=j+1;
if (r <= cell_fraction[i]) {cell_number_def[i]=cell_number_def[i]+1; break}
if (j==length(remaining_cells_shuffled)){cell_number_def[i]=cell_number_def[i]+1;}
}
}


cell_number=cell_number_def
print("Sum of cell numbers after bottleneck")
print(sum(cell_number))


end.time <- Sys.time()
time.taken <- end.time - start.time
print("Time to finish one week of experiment")
print(time.taken)

print(paste0("Number of extinct cell lineages by the end of week ", weeks))
print(sum(cell_number==0)) ## Number of extinct cell lineages
remaining_cell_lineages[weeks] <- sum(cell_number!=0)

cell.number.end.of.week <- cell_number
print("Correlation in cell numbers by cell identity between weeks")
print(cor.test(cell.number.end.of.week,cell.number.start.of.week))

}
remaining_cell_identity[[simulations]] <- cell_number
remaining_cells_count[simulations] <- sum(cell_number!=0)
#plot(remaining_cells_count)

remaining_cells_by_week_by_simulations <- rbind(remaining_cells_by_week_by_simulations, remaining_cell_lineages)

write.table(remaining_cells_by_week_by_simulations, gsub("xx", simulations, "$Address/remaining_cells_by_week_by_simulation_xx"))
write.table(cell_number, gsub("xx", simulations, "$Address/Analysis/cell_number_xx"))
### Writing the number of remaining cells and the frequency distribution stored in "cell number" ###

}

