## Atomistic simulation of protein evolution reveals sequence covariation and time-dependent fluctuations of site-specific substitution rates

This repository contains information about the computational procedures used to calculate the data for the manuscript

[Atomistic simulation of protein evolution reveals sequence covariation and time-dependent fluctuations of site-specific substitution rates](https://www.biorxiv.org/content/10.1101/2022.06.01.494278v1)

The repo focuses on aspects related to RosettaEvolve. RosettaEvolve is a set of methods used to simulated the evolution of proteins using an atomistic energy function and a stability-based fitness function. The code for RosettaEvolve is distrubuted through [Rosetta](https://www.rosettacommons.org/software) macromolecular modeling package, and together with the scripts found in this repo. 


There are four types of Rosetta evolve simulations presented in the study:

1) Sequence equilibration simuations used to equilibrate the sequence for a given energy offset value. Example data and scripts in the directory **equilibration_run**. 
2) Production runs simulating protein evolution along a single branch. The starting sequences and structures taken from the the end of the equilibration runs at given offsets. Example data and scripts in the directory p**roduction_run**.
3) Mutational trajectories where the trajectory is finished once a  single mutation is inserted. This is repeated as many times as there are sites in the protein. Example data and scripts in the directory **single_mutation_trajectory**. 
4) Simulation of a phylogenetic tree using RosettaEvolve. A tree crawling alogorithm recursively visits all the branches of a pre-defined phylogenetic tree. Example data and scripts in the directory **tree_crawling**.

Amino-acid site-rates are evalulated for some trajectories using a script based on a markov state model presented in the paper [A thermodynamic model of protein structure evolution explains empirical amino acid substitution matrices](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.4155)  The script is found in the directory **site_rates**.

More information can be found in README files in each of the subdirectories. 
