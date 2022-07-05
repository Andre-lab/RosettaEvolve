This repository contains information about the computational procedures used to calculate the data for the manuscript, focusing on the use of RosettaEvolve.

"Atomistic simulation of protein evolution reveals sequence covariation and time-dependent fluctuations of site-specific substitution rates", Christoffer Norn and Ingemar Andr√.

RosettaEvolve is a set of methods used to simulated the evolution of proteins using an atomistic energy function and a stability-based fitness function. There are four types of Rosetta evolve simulations 

1) Sequence equilibration simuations used to equilibrate the sequence for a given energy offset value. Example data and scripts in the directory equilibration_run. 
2) Production runs simulating protein evolutionary along a single branch. The starting sequences and structures taken from the the end of the equilibration runs at given offsets. Example data and scripts in the directory production_run.
3) Mutational trajectories where the trajectory is finished once a  single mutation is inserted. This is repeated as many times as there are sites in the protein. Example data and scripts in the directory single_mutation_trajectory. 
4) Simulation of a phylogenetic tree using RosettaEvolve. A tree crawling alogorithm recursively visits all the branches of a pre-defined phylogenetic tree. Example data and scripts in the directory tree_crawling.

Amino-acid site-rates are evalulated for some trajectories using a script based on a markov state model presented in "Norn C, Andre I, Theobald DL. A thermodynamic model of protein structure evolution explains empirical amino acid substitution matrices. Protein Sci. 2021;30(10):2057-68. Epub 2021/07/05. doi: 10.1002/pro.4155." The script is found in the directory site_rates.
