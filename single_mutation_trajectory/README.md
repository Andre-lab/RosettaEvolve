This directory contains an example xml script (production_run.xml) for evolving a sequence along a single branch. It does the following steps:

1) evolves the sequence for a certain number of trials
2) It dumps a PDB file
3) The effect of all point mutations (ddG) at all sites in the protein is calculated and stored in a text file.

1)-3) is carried out a certain number of times to allow for many multiple substutions occuring at the same site. In the manuscript on average 10 mutations per site is simulated.

The command to run one offset value is found in run_production.sh. The result of one run with an single offset is found in the directory example_ranks_files.

Comments on the xml script:

- steepness is the factor that convert the Rosetta energy (for the beta_nov16_cart energy function) onto the approximate kcal/mol scale (beta * T) in the Boltzmann function).  
- The NucleotideMutation move measures the ddG values at each position of the protein for all point mutations
- The EvolutionaryDynamics mover runs the evolutionary dynamics simulations and evolves the sequence
