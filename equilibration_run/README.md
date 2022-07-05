This directory contains an example xml script (equilibrate_sequence.xml) for equilibrating a PDB for further evolutionary simulations. It does the following steps:

1) evolves the sequence for a certain number of trials
2) It dumps a PDB file
3) The effect of all point mutations at all sites in the protein is calculated and stored in a text file. ddG can be calculated from this with the offset.

1)-3) is carried out a certain number of times to allow for many multiple substutions occuring at the same site. In the manuscript on average 10 mutations per site is simulated.

The command to run one offset value is found in run_equilibrate.sh. The result of one run with an single offset is found in the directory example_ranks_files.

Comments on the xml script:

- steepness is the factor that convert the Rosetta energy (for the beta_nov16_cart energy function) onto the approximate kcal/mol scale (beta * T) in the Boltzmann function).  
- The NucleotideMutation move measures the energy values at each position of the protein for all point mutations (from which the ddG can be calculated).
- The EvolutionaryDynamics mover runs the evolutionary dynamics simulations and evolves the sequence
