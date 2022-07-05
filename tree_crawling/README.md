A phylogenetic tree is simulated using a tree crawling simulation that simulates all branches along a pre-defined tree. We use a phylogenetic tree from a natural sequence alignment in this study.

The tree crawling algorithm recursively visits all nodes of a tree and has two parameter inputs.

python rosetta_evolution_wrapper_mpi.py --offset  -422.566  --mutationRate  4.07757e-06

-offset is the energy offset the evolutionary trajectory should be run at.
- mutationRate is the number of accepted mutations (non-synonymous) per trial. This is combined with the branch lenght stored in the phylogenetic tree to determine how many trials should
be attempted for each branch of the tree. Due to stochasticity, this will not give exactly the same branch length as the original tree. The mutationRate parameter is estimated from
the production run by evaluating the number of accepted mutations per trials in the trajectory progress file. See example in the production run example directory.

There are a couple of hard coded parameters within the script:
- the phylogenetic tree in Newick format
- The number of cores used in MPI
- Protein length and values used to give an estimation of the computation time used for the tree crawling.

Example of resulting trajectory data found in the directory trajectory.

The script calls and xml script that runs the RosettaEvolve simulation for each branch, simulate_branch_expected.xml, found here.
