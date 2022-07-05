#Example for running a production evolution run
rosetta_scripts.default.linuxgccrelease -s in/5azu.pdb -parser:protocol equilibrate_sequence.xml @flags -overwrite -out:prefix pdbs/5azu. -parser:script_vars id=5azu e_delta=-427.566 offset=-427.566 n_trials=127
