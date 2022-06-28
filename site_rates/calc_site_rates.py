import numpy as np
import mpmath
import matplotlib as mpl
import getpass
import glob

import sys
import time
import pickle
from optparse import OptionParser



##################################
# MISC
##################################
aas = "ARNDCQEGHILKMFPSTWYV*"
aas_nostop = 'ARNDCQEGHILKMFPSTWYV'
amino_acids_codon_no_stop = np.array([x for x in 'FFLLSSSSYYCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'])
amino_acids_codon = np.array([x for x in 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'])
stop_idx = np.where(amino_acids_codon == '*')[0]
aa2codonIdx_nostop = {aa: np.where(amino_acids_codon_no_stop == aa)[0] for aa in aas_nostop}
n_atoms = [89, 174,132,133,121,146, 147,75, 155, 131, 131,146,149,165,115,105,119,204,181,117]

aas_rev_no_stop = "VYWTSPFMKLIHGEQCDNRA"
codon_count_list = np.array([4, 6, 2, 2, 2, 2, 2, 4, 2, 3, 6, 2, 1, 2, 4, 6, 4, 1, 2, 4])

to_c = np.tile(codon_count_list, (20,1))
from_c = to_c.transpose()
relative_codon_rate = to_c
proposal_rates = codon_count_list/float(np.sum(codon_count_list))

aas_list = [x for x in aas_nostop]
aas_list_rev = [x for x in aas_rev_no_stop]
aa2no = {}
no2aa = {}
for i, aa in enumerate(aas):
    aa2no[aa] = i
    no2aa[i] = aa

# Setup index map for F61 to F20 conversion
F20_2_F61_idxs = np.empty((20, 20), dtype=object)
for j, J in enumerate(aas_nostop):
    for i, I in enumerate(aas_nostop):
        if j == i:  # We only considerkn the off diagonal elements...
            continue
        J_codons = np.where(amino_acids_codon_no_stop == J)[0]
        I_codons = np.where(amino_acids_codon_no_stop == I)[0]
        idxs_from = []
        idxs_to = []

        #print("correct list")
        for u in I_codons:
            for v in J_codons:
                idxs_from.append(v)
                idxs_to.append(u)
        idxs = (np.array(idxs_from), np.array(idxs_to))
        F20_2_F61_idxs[i][j] = idxs

np.set_printoptions(precision=3, edgeitems=64, linewidth=1000, threshold=1000, suppress=True)

# Vectorize the high precision exp function for numpy
def f(x):
    return mpmath.exp(x)
bigfloatEXP = np.vectorize(f)


########################################################################
# Data loading
########################################################################

def load_rank_file(filename, offset, score_style=None):
    data = {}
    chosen_energies_d = {}
    map_scanRes2pdbres = {}
    native_sequence = {}
    pop_these = []

    with open(filename,'r') as f_open:
        for line in f_open:
            ls = line.split()
            resi = int(ls[0])
            aa = ls[1]
            energy = float(ls[2])
            chosen = bool(int(ls[3])==1)
            resi0 = resi - 1

            if chosen:
                native_sequence[resi0] = aa
            if chosen and resi0 in chosen_energies_d and resi0 not in pop_these:
                pop_these.append(resi0)
            if chosen and resi0 not in chosen_energies_d:
                chosen_energies_d[resi0] = energy
            if resi0 not in data:
                data[resi0] = {}
            if aa not in data[resi0]:
                data[resi0][aa] = {}
            data[resi0][aa]['energy'] = energy


    #Find residues where we don't have all 20 amino acids
    for resi in data:
        if len(data[resi]) != 20:
            print("for", resi, "there are only", len(data[resi]), " amino acids with measured energy")
            for aa in "ARNDCQEGHILKMFPSTWYV":
                if aa not in data[resi]:
                    print("missing energies for ", resi+1, aa)

            if resi in pop_these:
                print("Already popping this one")
            else:
                pop_these.append(resi)

            sys.exit("Make sure everything is ok!!")


    # Setup map
    n_scanned = 0
    for n_pdb in range(0, len(data)):
        if n_pdb in pop_these:
            continue
        map_scanRes2pdbres[n_scanned] = n_pdb
        n_scanned += 1

    # Next remove those entries with multiple chosen amino acids (these are disulf that Rosetta could not handle)
    for resi in pop_these:
        if resi not in chosen_energies_d and resi in data:
            print("for", pdbid, resi, "we only popped res from the data list (means we could not measure energy of chosen aa)")
        if resi in data:
            data.pop(resi)
        if resi in chosen_energies_d:
            chosen_energies_d.pop(resi)
        if resi in native_sequence:
            native_sequence.pop(resi)

    # Place the data into 20x20 codon matrix instead
    M_e20 = np.zeros((len(data),20))
    native_sequence_out = []
    for i, resi in enumerate(data):
        # Because the protein is not relaxed completely different amount of energy
        # can be relax depending of which residue is being scanned. We account for this
        # by substracting each row in the matrix by the energy of the chosen amino
        # acid, so that each row will have 0.0, where the chosen amino acid is.
        if score_style=='nativeQinference':
            energies20 = np.array([data[resi][aa]['energy'] for aa in aas_nostop]) - chosen_energies_d[resi]
        else:
            energies20 = np.array([data[resi][aa]['energy'] for aa in aas_nostop]) - offset
        M_e20[i] = energies20
        native_sequence_out.append(native_sequence[resi])

    # Find the energy difference between the chosen aa and the best for that positions
    dist2bestAA = np.min(M_e20, axis=1)

    return chosen_energies_d, map_scanRes2pdbres, M_e20, dist2bestAA, native_sequence_out

def load_energy_data(data_folder='ranks/', offset=-300.0, score_style=None, rank_file=''):
#    f_list = glob.glob(data_folder + '/bl_*_*_*.txt')
    f_path = rank_file
    non_redundant_list = ["5azu"] ##

    f_list_d20 = {}
    chosen_energies = {}
    dist2bestAA = {}
    map_scanRes2pdbres = {}
    native_sequences = {}
    pdbid = rank_file.split('/')[1].split('_')[4].replace('.txt','' )

    if pdbid in non_redundant_list:
        chosen_energies[pdbid], map_scanRes2pdbres[pdbid], f_list_d20[pdbid], dist2bestAA[pdbid], native_sequences[pdbid] = load_rank_file(f_path, offset, score_style=score_style)

    M_e_all20 = np.zeros((0, 20))

    M_e_all_index_2_pdbid_and_res = {}
    pdbid_and_res_2_M_e_all_index = {}
    c = 0
    for pdbid in non_redundant_list:
        if pdbid not in pdbid_and_res_2_M_e_all_index:
            pdbid_and_res_2_M_e_all_index[pdbid] = []
        n = len(f_list_d20[pdbid])
        for resi, i in enumerate(range(c, c + n)):
            pdbid_and_res_2_M_e_all_index[pdbid].append((resi, i))
            M_e_all_index_2_pdbid_and_res[i] = (pdbid, map_scanRes2pdbres[pdbid][resi], resi)  # This map is 0 indexed in pdb
        c += n
        M_e_all20 = np.concatenate((M_e_all20, f_list_d20[pdbid]), axis=0)

    return M_e_all20


########################################################################
# Analysis functions
########################################################################

def Qnorm_from_Q(Q20, skip_list=None):
    if skip_list is not None:
        delete_list = []
        for aa in skip_list:
            delete_list.append(aa2no[aa])
        Q20 = np.delete(Q20, delete_list, axis=1)
        Q20 = np.delete(Q20, delete_list, axis=0)

    size = Q20.shape[0]
    # Set row sum
    d_index = np.diag_indices(size)
    Q20[d_index] = 0.0
    row_sums = np.nansum(Q20, axis=1)
    Q20[d_index] = -row_sums

    # Find pi and norm the matrix
    pi = np.ones(size).dot(np.linalg.inv(Q20 + np.ones((size, size))))
    relative_rate = -np.sum([Q20[i][i] * pi[i] for i in range(0, size)])
    Q20_normed = Q20 / relative_rate

    R_asym = Q20_normed / pi
    R_sym = (R_asym + R_asym.T) / 2
    R_sym[np.diag_indices(size)] = 0.0

    return Q20_normed, pi, R_sym

def get_pi20(w_matrix20, N=10000.0):
    exps = codon_count_list * bigfloatEXP(4 * N * w_matrix20)
    Zs = exps.sum(axis=1)
    pi20 = exps / Zs[:, None]
    pi20 = pi20.astype(float)
    return pi20

def p20_2_p61(pi20):
    pi20_codon_count_weighted = pi20 / codon_count_list
    pi61 = np.zeros((pi20.shape[0], 61))
    pi61[idxs_codon_vector_map] = pi20_codon_count_weighted[idxs_aa_vector_map]
    return pi61

def s2fixp(sel_coef_matrix, N=10000.0):
    fixp_M = (np.expm1(-2 * sel_coef_matrix)) / (np.expm1(-4 * N * sel_coef_matrix))

    # Correct to 1/2N where sel==0 (diagonal + elsewhere)
    indices_where_0 = np.where(sel_coef_matrix == 0)
    fixp_M[indices_where_0] = 1 / (2 * N)

    # Correct where nan (that would be when fixp is 0.)
    indices_where_nan = np.where(np.isnan(fixp_M))
    fixp_M[indices_where_nan] = 0.0

    return fixp_M

def hamm_dist(codon1, codon2):
    dist = 0
    for nt1, nt2 in zip(codon1, codon2):
        if nt1 != nt2:
            dist += 1
    return dist

def K80_rate(c1, c2, Q_nt_params):
    k = Q_nt_params['K80']
    nt2no = {'A': 0, 'G': 1, 'C': 2, 'T': 3}
    #             A  G  C  T
    Q = np.array([[0, k, 1, 1],
                 [k, 0, 1, 1],
                 [1, 1, 0 ,k],
                 [1, 1, k, 0]])
    product_rate = 1
    for nt1, nt2 in zip(c1, c2):
        if nt1 != nt2:
            product_rate *= Q[nt2no[nt1]][nt2no[nt2]]

    return product_rate

def make_Q_proposal(Q_nt_params, usePiCodonBias=False):
    bases = ['T', 'C', 'A', 'G']
    probability_of_mutation = Q_nt_params['p_codon_mut'] # 0.001
    codons = [a + b + c for a in bases for b in bases for c in bases]
    codon_connect = np.zeros((64, 64))

    for i, from_codon in enumerate(codons):
        for j, to_codon in enumerate(codons):
            codon_distance = hamm_dist(from_codon, to_codon)
            if codon_distance == 1:
                if 'K80' in Q_nt_params:
                    codon_connect[i][j] = probability_of_mutation + K80_rate(from_codon, to_codon, Q_nt_params)
            elif codon_distance == 2:
                 codon_connect[i][j] = probability_of_mutation #* K80_rate(from_codon, to_codon, Q_nt_params)
            elif codon_distance == 3:
                 codon_connect[i][j] = probability_of_mutation #* probability_of_mutation #* K80_rate(from_codon, to_codon, Q_nt_params)

    d_index = np.diag_indices(64)
    codon_connect[d_index] = 0.0
    row_sums = np.nansum(codon_connect, axis=1)
    codon_connect[d_index] = -row_sums

    codon_connect = np.delete(codon_connect, stop_idx, axis=1)
    codon_connect = np.delete(codon_connect, stop_idx, axis=0)

    return codon_connect

def get_F20(Q61, pi_codons):
    F20 = np.zeros((20,20))
    F61 = pi_codons[:, np.newaxis] * Q61
    for j, J in enumerate(aas_nostop):
        for i, I in enumerate(aas_nostop):
            if j == i:  # We only consider the off diagonal elements...
                continue
            F20[i][j] = F61[F20_2_F61_idxs[i][j]].sum()

    site_flux = F20.sum()

    d_index = np.diag_indices(20)
    F20[d_index] = 0.0
    row_sums = np.sum(F20, axis=1)
    F20[d_index] = -row_sums

    return F20, site_flux


def convertE2substitutionMatrix(EM20, Q_nt_params, s=1.68, offset=-300, N=10000.0, usePiCodonBias=False, protein_specific_Ns=None, protein_specific_Os=None, NstDev=0.0, score_style=None):
    protein_length = EM20.shape[0]

    proposal_Q = make_Q_proposal(Q_nt_params, usePiCodonBias=usePiCodonBias) # This includes K80 rates with k=2.0.

    if score_style=='nativeQinference':
        w_matrix20 = 1/(1 + np.exp((EM20 - offset) * s)) # n_resi, 64
    else:
        w_matrix20 = 1/(1 + np.exp(EM20 * s)) # n_resi, 64

    w_to = np.array([np.tile(w_vector, (20,1)) for w_vector in w_matrix20])
    w_from = np.array([m.transpose() for m in w_to])
    sel_coef_matrix = (w_to - w_from)/w_from

    fixp_M = s2fixp(sel_coef_matrix, N=N)
    pi20 = get_pi20(w_matrix20, N=N)

    # Remap to a 61x61 codon matrix
    fixp_M61 = np.zeros((protein_length, 61, 61))
    fixp_M61[idxs_codon_matrix_map] = fixp_M[idxs_aa_matrix_map]
    pi61 = p20_2_p61(pi20)

    # Multiply on the proposal frequencies
    codon_Q_M = proposal_Q * fixp_M61

    # Set diagonal elements to 0
    [np.fill_diagonal(M, 0) for M in codon_Q_M]
    rates_nt = [np.sum(pi61[i][:, np.newaxis] * codon_Q_M[i] ) for i in range(0, len(codon_Q_M))]
    # Set diagonal elements to -row sum
    [np.fill_diagonal(M, -np.nansum(M, axis=1)) for M in codon_Q_M]

    # Convert the Q64 matrices to Q20 and find the aa stationary freqiencies
    Fs_return = []
    rates = []
    n = 1
    for i in range(0, protein_length):
        n += 1
        F, rate_aa = get_F20(codon_Q_M[i], pi61[i])
        rates.append(rate_aa)
        Fs_return.append(F)

    return np.array(Fs_return), pi20, np.array(rates), np.array(rates_nt), codon_Q_M

def Q_norm(Q, pi):
    size = Q.shape[0]

    d_index = np.diag_indices(size)
    Q[d_index] = 0.0
    row_sums = np.nansum(Q, axis=1)
    Q[d_index] = -row_sums

    diagonals = np.diagonal(Q)
    rate = -np.sum(diagonals * pi)

    return Q / rate

def flux2Q(M_pi, Fs, rates, outstr, cutoff=1e-300):
    pi_mean_unfiltered = np.nansum(M_pi, axis=0)
    pi_mean_unfiltered = pi_mean_unfiltered / np.nansum(pi_mean_unfiltered)
    nans = []
    n_rates = 0
    n_no_rate = 0
    # D2Vs = []
    pi_mean = np.zeros(20)
    F_distribution = np.zeros((len(Fs), 20, 20))
    for i, F in enumerate(Fs):
        if rates[i] > cutoff:
            n_rates += 1
            F_normed = F / rates[i]
            F_distribution[i, :, :] = F_normed
            pi_mean += M_pi[i]
        else:
            n_no_rate += 1
            nans.append(i)
    pi_mean = pi_mean / np.sum(pi_mean)

    Q_mean = np.nansum(F_distribution, axis=0) / pi_mean[:, np.newaxis]
    Q_mean_normed = Q_norm(Q_mean, pi_mean)

    return pi_mean_unfiltered, pi_mean, Q_mean_normed, F_distribution, n_no_rate, n_rates

def precompute_61to20_map():
    # Create map over (resi, 20 long w/e) to (resi, 61 long w/e)
    # noinspection PyUnboundLocalVariable
    try:
        #print("loaded pickled maps")
        idxs_aa_matrix_map = pickle.load(open('idxs_aa_matrix_map.pickle', 'rb'))
        idxs_codon_matrix_map = pickle.load(open('idxs_codon_matrix_map.pickle', 'rb'))
        idxs_aa_vector_map = pickle.load(open('idxs_aa_vector_map.pickle', 'rb'))
        idxs_codon_vector_map = pickle.load(open('idxs_codon_vector_map.pickle', 'rb'))
    except:
        t0 = time.time()
        # Make matrix maps
        resi_idx = []
        from_idx_61 = []
        to_idx_61 = []
        from_idx_20 = []
        to_idx_20 = []
        for resi in range(0, M_e_all20.shape[0]):
            for from_idx_aa, from_aa in enumerate(aas_nostop):
                from_idxs_codon = aa2codonIdx_nostop[from_aa]
                for to_idx_aa, to_aa in enumerate(aas_nostop):
                    to_idxs_codon = aa2codonIdx_nostop[to_aa]
                    for from_c_idx in from_idxs_codon:
                        for to_c_idx in to_idxs_codon:
                            resi_idx.append(resi)
                            from_idx_61.append(from_c_idx)
                            to_idx_61.append(to_c_idx)
                            from_idx_20.append(from_idx_aa)
                            to_idx_20.append(to_idx_aa)
        idxs_aa_matrix_map = (np.array(resi_idx), np.array(from_idx_20), np.array(to_idx_20))
        idxs_codon_matrix_map = (np.array(resi_idx), np.array(from_idx_61), np.array(to_idx_61))

        # Make vector maps
        resi_idxs = []
        idxs_61 = []
        idxs_20 = []
        for resi in range(0, M_e_all20.shape[0]):
            for aa_idx, from_aa in enumerate(aas_nostop):
                idxs_codon = aa2codonIdx_nostop[from_aa]
                for c_idx in idxs_codon:
                    idxs_20.append(aa_idx)
                    idxs_61.append(c_idx)
                    resi_idxs.append(resi)
        idxs_aa_vector_map = (np.array(resi_idxs), np.array(idxs_20))
        idxs_codon_vector_map = (np.array(resi_idxs), np.array(idxs_61))

        with open('idxs_aa_matrix_map.pickle', 'wb') as f:
            pickle.dump(idxs_aa_matrix_map, f)
        with open('idxs_codon_matrix_map.pickle', 'wb') as f:
            pickle.dump(idxs_codon_matrix_map, f)
        with open('idxs_aa_vector_map.pickle', 'wb') as f:
            pickle.dump(idxs_aa_vector_map, f)
        with open('idxs_codon_vector_map.pickle', 'wb') as f:
            pickle.dump(idxs_codon_vector_map, f)

        t1 = time.time()
        print("made maps")
        print(t1 - t0, "seconds")

    return idxs_aa_matrix_map, idxs_codon_matrix_map, idxs_aa_vector_map, idxs_codon_vector_map


################################
# Main
################################

### Build optionparser ###
parser = OptionParser(usage="usage: %prog [options] FILE", version="0.1")
parser.add_option("-r", "--rank_file", type="str", dest="rank_file", metavar="STR", help="rank file")
parser.add_option("-o", "--offset", type="float", dest="offset", metavar="FLOAT", help="offset")
#parser.add_option("-s", "--score_style", type="str", dest="score_style", metavar="STR", help="score styl")
(opts, args) = parser.parse_args()
parser.set_defaults(rank_file='')

#offset = 12.5 # XXX
#offset = -332.870 +12.5 # XXX
offset = opts.offset
Q_nt_params = {'K80': 2.7, 'p_codon_mut': 0.1} # XXX
#N = 10**2.2 # XXX
N = 15848
s = 0.86 # XXX
outstr = 'test_for_dump_1' # XXX
#score_style = 'nativeQinference'
score_style = ''

M_e_all20 = load_energy_data(offset=offset, score_style = score_style, rank_file=opts.rank_file)
#print(M_e_all20)
idxs_aa_matrix_map, idxs_codon_matrix_map, idxs_aa_vector_map, idxs_codon_vector_map = precompute_61to20_map()
Fs, M_pi, rates, rates_nt, codon_Q_M = convertE2substitutionMatrix(M_e_all20, Q_nt_params, offset=offset, s=s, N=N, score_style=score_style)
pi_mean_unfiltered, pi_mean, Q_mean_normed, F_distribution, n_no_rate, n_rates = flux2Q(M_pi, Fs, rates, outstr)

f = opts.rank_file
f_rank = f.split('/')[1]
base = f_rank.replace('.txt','')
f_out = open('rates/'+base+'.rates',"w")

# This is rates
for i,r in enumerate(rates):
#    print(i, r*N)
       l = str(i) + ' ' + str(r*N) + '\n'
       f_out.write(l)

# This is the Q matrix
#print(Q_mean_normed)

