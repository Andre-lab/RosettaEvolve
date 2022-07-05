Script for calculating site mutation rates from calculated ddG values at each site of a protein

The script calc_site_rates.py takes two inputs:

- A rank file that stores all the ddG values
- A values of the offset used to simulate the protein. The offset is used to calculate ddG values from the absolute energy measurments in the rank file

In addition to this there are a couple of hard coded values in the script, kappa, rho, effective population size (Neff) and the factor that converts the Rosetta energies into kcal/mol scale (s).
These values can be modified in the script.
