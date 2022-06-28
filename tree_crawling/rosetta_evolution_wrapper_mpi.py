__author__ = 'christoffernorn'

from Bio import Phylo
from cStringIO import StringIO
import networkx #, pylab
import dendropy
import subprocess as subprocess
import shutil
import os
import time
import sys
# tree = Phylo.read('example.dnd', 'newick')
from optparse import OptionParser

## Test that correct versions are being
if sys.version_info[0] != 2:
	raise Exception("Must be using Python 2.x")
if networkx.__version__ != '1.11':
	raise Exception("Please use network 1.11")


### Build optionparser ###
parser = OptionParser(usage="usage: %prog [options] FILE", version="0.1")
parser.add_option("-r", "--mutationRate", type="float", dest="mutationRate", metavar="FLOAT", help="Mutation rate")
parser.add_option("-o", "--offset", type="float", dest="offset", metavar="FLOAT", help="energy offset")

(opts, args) = parser.parse_args()
parser.set_defaults(mutationRate=1.0, offset=0.0)

# clean up
shutil.rmtree('traj')
shutil.rmtree('pdbs')
shutil.rmtree('logs')
# Setup
os.makedirs('traj')
os.makedirs('pdbs')
os.makedirs('logs')

def recurse_tree(G, from_node, to_node, is_dry_run):
    n_mutations_to_introduce = str(G.get_edge_data(to_node, from_node)['weight'] * protein_length )

    if is_dry_run:
        cmd_string = "sleep " + str((float(n_mutations_to_introduce)/100))
    else:
        pdb_in = "pdbs/" + str(from_node) + ".pdb"
        pdb_out = "pdbs/" + str(to_node) + "_"
        trajectory_file = "traj/" + str(from_node) + "_" + str(to_node) + ".traj"
        exe = "/pfs/nobackup/home/i/iandre/rosetta_evolve/Rosetta/rosetta_scripts.static.linuxgccrelease "
        cmd_string = exe + '-s ' + pdb_in + ' @flags ' + ' -out:prefix ' + pdb_out \
                    + ' -parser:script_vars progress_id=' + trajectory_file \
                    + ' -parser:script_vars branch_length=' + n_mutations_to_introduce \
                    + ' -parser:script_vars offset=' + str(opts.offset) \
                    + ' -parser:script_vars mutationRate=' + str(opts.mutationRate)
        logs=True
        if logs:
            cmd_string += " > logs/"+str(from_node)+"_"+str(to_node)+".log"

        global mutations_so_far
        mutations_so_far += int(round((G.get_edge_data(to_node, from_node)['weight'] * protein_length)))
        print str(round(float(mutations_so_far)/float(total_number_of_mutations)*100,1))+"% complete."+" Simulating trajectory between " + str(from_node) + " and " + str(to_node) + " with distance is " + str(
                G.get_edge_data(to_node, from_node)['weight']) + " with branch length of " + n_mutations_to_introduce + " follow trajectory at " + trajectory_file

    print "Opening process: "+cmd_string
    undergoing_work_paths[(from_node,to_node)] = subprocess.Popen(cmd_string, shell=True)

    wait = True
    while wait:
        # Let's check whether any trajectories are done and then update the paths_ready_to_explore accordingly
        if len(undergoing_work_paths) > 0:
            for path in list(undergoing_work_paths):
                is_done = (undergoing_work_paths[path].poll() == 0)
                if is_done:
                    undergoing_work_paths.pop(path)
                    from_node = path[0]
                    to_node = path[1]
                    known_node_list.append(to_node)
                    newly_discovered_path_frontier = [(to_node,x) for x in net.neighbors(to_node) if x not in known_node_list]
                    for (x, y) in newly_discovered_path_frontier:
                        paths_ready_to_explore.append((x,y))

                    # and clean up any output files
                    pdb_out_file_name = str(to_node) + "_" + str(from_node) + "_0001.pdb"
                    if is_dry_run:
                        continue
                    shutil.move("pdbs/" + pdb_out_file_name, "pdbs/" + str(to_node) + ".pdb")

        # If we have available CPUs and there is something to explore, we do that
        if len(paths_ready_to_explore) > 0 and len(undergoing_work_paths) < max_CPU:
            from_node = paths_ready_to_explore[0][0]
            to_node = paths_ready_to_explore[0][1]
            paths_ready_to_explore.pop(0)
            recurse_tree(G, from_node, to_node, is_dry_run)

        # If there is no current jobs and nothing ready to explore, then we must be done
        if len(paths_ready_to_explore) == 0 and len(undergoing_work_paths) == 0:
            wait = False

        # wait a bit
        time.sleep(0.001)


#############################################################
# Do for a test tree
#############################################################

#
# max_CPU = 1
# treedata = "(A:3,B:5,(C:6,D:1)E:2)F"
# handle = StringIO(treedata)
# tree = Phylo.read(handle, "newick")
# net = Phylo.to_networkx(tree)
#
# from_node = tree.find_clades("E").next()
# to_node = net.neighbors(from_node)[0]
#
# paths_ready_to_explore = [(from_node, to_n) for to_n in net.neighbors(from_node)[1::]]
# known_node_list = [from_node]
# undergoing_work_paths = {}
#
# recurse_tree(net, from_node, to_node)


# ############################################################
# Commandline args
# ############################################################

treefile = "in/tree.newick"
pdbfile = "in/in.pdb"
runid = str(1)
max_CPU = 28
protein_length = 127

time_per_mutation = 181 # seconds

#############################################################
# Do for a real tree
#############################################################
# First relabel the nodes
filepath = treefile  # ../build_trees/RAxML_bestTree.PF00381.align_raxml
tree = dendropy.Tree.get(path=filepath, schema="newick")
treelen = tree.length()
leaves = tree.leaf_nodes()
leaf_num = len(leaves)
time_per_branch_length = 1
print "The total branch length is: %f" % treelen
print "Ideally this will take " + str((treelen * protein_length * time_per_mutation)/float(max_CPU)) + " seconds"
print "Will generate %d extant sequences" % leaf_num

tree.deroot()
tree.update_bipartitions()

# write node number to node label
i = leaf_num
for node in tree.preorder_node_iter():
    if node.is_internal():
        node.label = "node_" + str(i)
        i += 1

# Convert the tree to phylo format
handle = StringIO(str(tree))
tree = Phylo.read(handle, "newick")
net = Phylo.to_networkx(tree)

# We might run into recursion depth problems if there are too many edges
graph_edges = len(net.edges())
if graph_edges > 1000:
    sys.setrecursionlimit(graph_edges * 2 )

print "I will be generated", graph_edges, "trajectories."

# Recurse the tree - dry mode (to check for potential issues)

skip_dry_check = True
if not skip_dry_check:
    center_node = networkx.center(net)[0]  # To optimize run time, we start from a node in the center of the graph
    from_node = tree.find_clades(center_node).next()
    shutil.copyfile(pdbfile,'pdbs/'+str(from_node)+'.pdb')
    to_node = net.neighbors(from_node)[0]
    paths_ready_to_explore = [(from_node, to_n) for to_n in net.neighbors(from_node)[1::]]
    known_node_list = [from_node]
    undergoing_work_paths = {}
    is_dry_run = True
    recurse_tree(net, from_node, to_node, is_dry_run)
    print "I checked that the tree can be recursed, and there are no problems..."
    print "If the CPU estimate looks reasonable change skip_dry_check to True and start the production run!"
    sys.exit()

# Recurse the tree -- production mode

#######################################################
# Add %complete counter
#######################################################

print "starting production run..."
start = time.time()
total_number_of_mutations = treelen * protein_length
mutations_so_far = 0
center_node = networkx.center(net)[0]  # To optimize run time, we start from a node in the center of the graph
from_node = tree.find_clades(center_node).next()
shutil.copyfile(pddbfile,'pdbs/'+str(from_node)+'.pdb')
to_node = net.neighbors(from_node)[0]
paths_ready_to_explore = [(from_node, to_n) for to_n in net.neighbors(from_node)[1::]]
known_node_list = [from_node]
undergoing_work_paths = {}
is_dry_run = False
recurse_tree(net, from_node, to_node, is_dry_run)
end = time.time()
print "Total run time ", end - start
