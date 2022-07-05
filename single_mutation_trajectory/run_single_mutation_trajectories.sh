#Example for running the single mutation run 
rosetta_scripts.default.linuxgccrelease -s in/5azu.pdb -parser:protocol single_mutation_only.xml @flags -overwrite -out:prefix pdbs/5azu. -parser:script_vars id=5azu e_delta=-427.566 offset=-427.566 n_trials=1
