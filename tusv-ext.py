#   author: Jesse Eaton, Xuecong Fu
#   the file is originated from tusv.py from TUSV by Jesse. Xuecong Fu fixed bugs and extend to current model TUSV-est.

# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #
import copy
import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import random
import numpy as np
import multiprocessing as mp

from datetime import datetime
from graphviz import Digraph
from ete2 import Tree          # for creating phylogenetic trees for .xml output
from Bio import Phylo          # for creating phylogenies to export as phylo .xml files
from cStringIO import StringIO # for converting string to file (for creating initial phylo .xml)

sys.path.insert(0, 'model/')
sys.path.insert(0, 'help/')
import solver as sv
import file_manager as fm      # sanitizes file and directory arguments
import generate_matrices as gm # gets F, Q, G, A, H from .vcf files
import printer as pt
import vcf_help as vh
import pickle
from snv_matching import snv_assign

# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

MAX_NUM_LEAVES = 10
MAX_COPY_NUM = 20
MAX_CORD_DESC_ITERS = 1000
MAX_RESTART_ITERS = 1000
NUM_CORES = mp.cpu_count()
METADATA_FNAME = 'data/2017_09_18_metadata.vcf'
STR_DTYPE = 'S50'


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
    args = get_args(argv)
    write_readme(args['output_directory'], args)
    unmix(args['input_directory'], args['output_directory'], args['num_leaves'], args['c_max'], args['lambda1'], args['lambda2'], args['restart_iters'], args['cord_desc_iters'], args['processors'], args['time_limit'], args['metadata_file'], args['num_subsamples'], args['overide_lambdas'], args['constant'], args['sv_upperbound'], args['only_leaf'], args['collapse'], args['threshold'], args['multi_num_clones'])


#  input: num_seg_subsamples (int or None) number of segments to include in deconvolution. these are
#           in addition to any segments contining an SV as thos are manditory for the SV. None is all segments
def unmix(in_dir, out_dir, n, c_max, lamb1, lamb2, num_restarts, num_cd_iters, num_processors, time_limit, metadata_fname, \
          num_seg_subsamples, should_overide_lambdas, const, sv_ub, only_leaf, collapse, threshold, multi_num_clones=False):
    print("unmix")

    F_phasing_full, F_unsampled_phasing_full, Q_full, Q_unsampled_full, G, A, H, bp_attr, cv_attr, F_info_phasing, \
    F_unsampled_info_phasing, sampled_snv_list_sort, unsampled_snv_list_sort, sampled_sv_list_sort, unsampled_sv_list_sort = gm.get_mats(in_dir, n, const=const, sv_ub=sv_ub)
    Q_full, Q_unsampled_full, G, A, H, F_phasing_full, F_unsampled_phasing_full = check_valid_input(Q_full, Q_unsampled_full,G, A, H, F_phasing_full, F_unsampled_phasing_full)

    np.savetxt(out_dir + "/F_info_phasing.csv", F_info_phasing, delimiter='\t', fmt='%s')
    np.savetxt(out_dir + "/F_unsampled_info_phasing.csv", F_unsampled_info_phasing, delimiter='\t', fmt='%s')
    np.savetxt(out_dir + "/sampled_snv_list_sort.csv", sampled_snv_list_sort, delimiter='\t', fmt='%d')
    np.savetxt(out_dir + "/unsampled_snv_list_sort.csv", unsampled_snv_list_sort, delimiter='\t', fmt='%d')
    np.savetxt(out_dir + "/sampled_sv_list_sort.csv", sampled_sv_list_sort, delimiter='\t', fmt='%d')
    np.savetxt(out_dir + "/unsampled_sv_list_sort.csv", unsampled_sv_list_sort, delimiter='\t', fmt='%d')
    F_phasing, Q, Q_unsampled, org_indxs = randomly_remove_segments(F_phasing_full, Q_full, Q_unsampled_full, num_seg_subsamples)
    np.savetxt(out_dir + '/F_phasing.tsv', F_phasing, delimiter='\t', fmt='%.8f')
    np.savetxt(out_dir + '/F_unsampled_phasing_full.tsv', F_unsampled_phasing_full, delimiter='\t', fmt='%.8f')
    # replace lambda1 and lambda2 with input derived values if should_orveride_lamdas was specified
    m = len(F_phasing)
    l_g, r = Q.shape
    g_un = Q_unsampled.shape[0]
    print('The num of features of F is '+str(l_g)+ ', the num of copy numbers is ' +str(r)+ ', the num of unsampled SNV is ' + str(g_un)+ '.')
    if should_overide_lambdas:

        lamb1 = float(l_g + 2*r) / float(2*r) * float(m) / float(2 * (n-1) )/2
        lamb2 = float(l_g + 2*r) / float(l_g)/2

    Us, Cs, Es, As, obj_vals, Rs, Ws, W_SVs, W_SNVs = [], [], [], [], [], [], [], [], []
    num_complete = 0
    if not multi_num_clones:
        for i in xrange(0, num_restarts):
            U, C, E, A_, R, W, W_SV, W_SNV, obj_val, err_msg = sv.get_UCE(F_phasing, Q, G, A, H, n, c_max, lamb1, lamb2, num_cd_iters, time_limit, only_leaf)
            printnow(str(i + 1) + ' of ' + str(num_restarts) + ' random restarts complete\n')
            Us.append(U)
            Cs.append(C)
            Es.append(E)
            As.append(A_)
            Rs.append(R)
            Ws.append(W)
            W_SVs.append(W_SV)
            W_SNVs.append(W_SNV)
            obj_vals.append(obj_val)

        best_i = 0
        best_obj_val = obj_vals[best_i]
        for i, obj_val in enumerate(obj_vals):
            if obj_val < best_obj_val:
                best_obj_val = obj_val
                best_i = i

        with open(out_dir + "/training_objective", 'w') as f:
            f.write(str(best_obj_val))
        E_pre = copy.deepcopy(Es[best_i])
        R_pre = copy.deepcopy(Rs[best_i])
        W_pre = copy.deepcopy(Ws[best_i])
        if collapse:
            U_best, C_best, E_best, A_best, R_best, W_best, W_SV_best, W_SNV_best = collapse_nodes(Us[best_i], Cs[best_i], Es[best_i], As[best_i], Rs[best_i], Ws[best_i], W_SVs[best_i], W_SNVs[best_i], threshold,only_leaf)
        else:
            U_best, C_best, E_best, A_best, R_best, W_best, W_SV_best, W_SNV_best = Us[best_i], Cs[best_i], Es[best_i], As[best_i], Rs[best_i], Ws[best_i], W_SVs[best_i], W_SNVs[best_i]
        min_node, min_dist, W_unsampled = snv_assign(C_best[:, -2*r:], Q_unsampled, A_best, E_best, U_best, F_unsampled_phasing_full)
        np.savetxt(out_dir + "/unsampled_assignment.csv", min_node, delimiter=',')
        np.savetxt(out_dir + "/unsampled_assignment_dist.csv", min_dist, delimiter=',')
        ### concatenate unsampled SV and SNV list
        W_SV_unsampled = W_unsampled[:,:len(unsampled_sv_list_sort)]
        W_SNV_unsampled = W_unsampled[:,len(unsampled_sv_list_sort):]
        W_con = concatenate_W(W_SV_best, W_SV_unsampled, W_SNV_best, W_SNV_unsampled, sampled_sv_list_sort, unsampled_sv_list_sort, sampled_snv_list_sort, unsampled_snv_list_sort)
        writer = None #build_vcf_writer(F_phasing_full, C_best, org_indxs, G, Q, bp_attr, cv_attr, metadata_fname)
        B = create_binary_matrix(W_con, A_best)
        write_to_files(out_dir, l_g, U_best, C_best, E_best, R_best, W_best, W_SV_best, W_SNV_best, W_unsampled, W_con, obj_vals[best_i], F_phasing_full, F_unsampled_phasing_full, org_indxs, writer, E_pre, R_pre, W_pre, B, A_best)
    else:
        training_obj = np.zeros(n-1)
        for n_ in range(2, n+1):
            U, C, E, A_, R, W, W_SV, W_SNV, obj_val, err_msg = sv.get_UCE(F_phasing, Q, G, A, H, n_, c_max, lamb1,
                                                                              lamb2, num_cd_iters, time_limit, only_leaf)
            printnow(str(n_) + ' of ' + str(num_restarts) + ' num of clones restarts complete\n')
            training_obj[n_-2] = obj_val
            E_pre = copy.deepcopy(E)
            R_pre = copy.deepcopy(R)
            W_pre = copy.deepcopy(W)
            if collapse:
                U, C, E, A_, R, W, W_SV, W_SNV = collapse_nodes(U,C,E,A_,R,W,W_SV, W_SNV,threshold,only_leaf)

            min_node, min_dist, W_SNV_unsampled = snv_assign(C[:, -2 * r:], Q_unsampled, A_, E, U,F_unsampled_phasing_full)
            np.savetxt(out_dir + "/unsampled_SNV_assignment.csv", min_node, delimiter=',')
            np.savetxt(out_dir + "/unsampled_SNV_assignment_dist.csv", min_dist, delimiter=',')
            W_con, W_snv_con = concatenate_W(W_SV, W_SNV, W_SNV_unsampled, sampled_snv_list_sort,
                                             unsampled_snv_list_sort)
            writer = build_vcf_writer(F_phasing_full, C, org_indxs, G, Q, bp_attr, cv_attr, metadata_fname)
            B = create_binary_matrix(W_con, A)
            if not os.path.exists(out_dir + '/num_clone_' + str(n_)):
                os.mkdir(out_dir + '/num_clone_' + str(n_))
            write_to_files(out_dir + '/num_clone_' + str(n_) + '/', l_g, U, C, E, R, W, W_SV, W_SNV, W_SNV_unsampled,W_con, obj_val, F_phasing_full,
                           F_unsampled_phasing_full, org_indxs, writer, E_pre, R_pre, W_pre, B, A_)
        np.savetxt(out_dir + '/training_obj_list.csv', training_obj, delimiter='\t')

def create_binary_matrix(W_con, A):
    B = copy.deepcopy(W_con)
    ad_pairs = np.where(A == 1)
    for i in range(len(ad_pairs[0])):
        ancestor_idx = ad_pairs[0][i]
        descendant_idx = ad_pairs[1][i]
        B[descendant_idx] += B[ancestor_idx]
    B[B > 1] = 1
    return B

# concatenating W matrix for SVs and SNVs
def concatenate_W(W_SV_TUSV, W_SV_MATCHING, W_SNV_TUSV, W_SNV_MATCHING, sampled_sv_list_sort, unsampled_sv_list_sort, sampled_snv_list_sort, unsampled_snv_list_sort):
    n, l_sampled = W_SV_TUSV.shape
    l_unsampled = W_SV_MATCHING.shape[1]
    l = l_sampled + l_unsampled
    g_sampled = W_SNV_TUSV.shape[1]
    g_unsampled = W_SNV_MATCHING.shape[1]
    g = g_sampled + g_unsampled
    W_con = np.zeros((n, l + g))
    W_snv_con = np.zeros((n, g))
    if l_unsampled != 0:
        W_con[:, sampled_sv_list_sort] = W_SV_TUSV
        W_con[:, unsampled_sv_list_sort] = W_SV_MATCHING
    else:
        W_con[:, :l]= W_SV_TUSV
    if g_sampled != 0:
        W_snv_con[:, sampled_snv_list_sort] = W_SNV_TUSV
        W_snv_con[:, unsampled_snv_list_sort] = W_SNV_MATCHING
    else:
        W_snv_con = W_SNV_MATCHING
    W_con[:, l:] = W_snv_con
    return W_con

# create tree from W matrix
def W2tree(W_sv_total, W_snv_total, E):
    edge_list = np.where(E == 1)
    tree = {}
    mutations = {}
    for i in range(len(edge_list[0])):
        parent = edge_list[0][i]
        child = edge_list[1][i]
        if parent not in tree.keys():
            tree[parent] = []
        tree[parent].append(child)
        mut_sv_list = list(np.where(W_sv_total[child, :] == 1))
        mutations[child] = ['sv_' + str(j) for j in mut_sv_list]
        mut_snv_list = list(np.where(W_snv_total[child, :] == 1))
        mutations[child] += ['snv_' + str(j) for j in mut_snv_list]
    return tree, mutations

# collapse nodes
def collapse_nodes(U, C, E, A, R, W, W_SV, W_SNV, threshold=0.0, only_leaf=False):
    print("Loading collapse nodes")
    # generate the tree
    tree = ModifyTree(E)
    if not only_leaf:
        # collapse the branches with 0 length
        branch_remove_idx = []
        for i in xrange(tree.N-1, -1, -1):
            for j in xrange(tree.N-1, -1, -1):
                if int(E[i, j]) == 1 and sum(W[j, :]) == 0 and R[i,j] == 0:
                    branch_remove_idx.append(j)
        for node in branch_remove_idx:
            target = tree.cp_tree[node]
            U[:, target] += U[:, node]
            if not tree.is_leaf(node):
                for child in tree.tree[node]:
                    R[target, child] = R[node, child]
            tree.delete_node(node)

        # collapse the nodes with 0 frequency
        freq_remove_idx = []
        freq_leaf_remove_idx = []
        for i in xrange(tree.N-1, -1, -1):
            if i in branch_remove_idx:
                continue
            if np.mean(U[:, i]) <= threshold:
                if tree.num_children(i) == 1:
                    freq_remove_idx.append(i)
                elif tree.is_leaf(i):
                    freq_leaf_remove_idx.append(i)
        for node in freq_remove_idx:
            target = tree.tree[node][0]
            parent = tree.cp_tree[node]
            tree.delete_node(node)
            W[target, :] += W[node, :]
            W_SV[target,:] += W_SV[node,:]
            W_SNV[target, :] += W_SNV[node, :]
            R[parent, target] = R[parent, node] + R[node, target]
        for node in freq_leaf_remove_idx:
            tree.delete_node(node)
    else:
        # collapse the branches with 0 length and the child of the branch doesn't belong to leaf nodes
        branch_remove_idx = []
        for i in xrange(tree.N - 1, -1, -1):
            for j in xrange(tree.N - 1, -1, -1):
                if int(E[i, j]) == 1 and sum(W[j, :]) == 0 and R[i, j] == 0 and not tree.is_leaf(j):
                    branch_remove_idx.append(j)
        for node in branch_remove_idx:
            target = tree.cp_tree[node]
            U[:, target] += U[:, node]
            if not tree.is_leaf(node):
                for child in tree.tree[node]:
                    R[target, child] = R[node, child]
            tree.delete_node(node)

        # collapse the leaf nodes with 0 frequency
        freq_remove_idx = []
        freq_leaf_remove_idx = []
        for i in xrange(tree.N - 1, -1, -1):
            if i in branch_remove_idx:
                continue
            if np.mean(U[:, i]) <= threshold and tree.is_leaf(i):
                freq_leaf_remove_idx.append(i)
        for node in freq_leaf_remove_idx:
            tree.delete_node(node)
        for i in xrange(tree.N - 1, -1, -1):
            if tree.num_children(i) == 1:
                freq_remove_idx.append(i)
        for node in freq_remove_idx:
            target = tree.tree[node][0]
            parent = tree.cp_tree[node]
            tree.delete_node(node)
            W[target, :] += W[node, :]
            W_SV[target, :] += W_SV[node, :]
            W_SNV[target, :] += W_SNV[node, :]
            R[parent, target] = R[parent, node] + R[node, target]


    # delete those nodes
    remove_idx = branch_remove_idx + freq_remove_idx + freq_leaf_remove_idx
    print('Nodes ', remove_idx, 'will be collapsed.')
    U_new = np.delete(U, remove_idx, axis=1)
    C_new = np.delete(C, remove_idx, axis=0)
    A_new = np.delete(A, remove_idx, axis=0)
    A_new = np.delete(A_new, remove_idx, axis=1)
    E_new = np.delete(tree.E, remove_idx, axis=0)
    E_new = np.delete(E_new, remove_idx, axis=1)
    R_new = np.delete(R, remove_idx, axis=0)
    R_new = np.delete(R_new, remove_idx, axis=1)
    W_new = np.delete(W, remove_idx, axis=0)
    W_SV_new = np.delete(W_SV, remove_idx, axis=0)
    W_SNV_new = np.delete(W_SNV, remove_idx, axis=0)
    print("collapse", U_new.shape, C_new.shape)
    return U_new, C_new, E_new, A_new, R_new, W_new, W_SV_new, W_SNV_new


class ModifyTree:
    def __init__(self, E):
        self.cp_tree = {}
        self.tree = {}
        self.E = E
        self.N = len(E)
        for i in xrange(self.N - 1, -1, -1):
            for j in xrange(self.N - 1, -1, -1):
                if int(E[i, j]) == 1:
                    self.cp_tree[j] = i
                    if i not in self.tree.keys():
                        self.tree[i] = [j]
                    else:
                        self.tree[i].append(j)

    def delete_node(self, idx):
        if self.is_root(idx):
            if self.num_children(idx) > 1:
                raise('Cannot delete root node with more than one child!')
            child = self.tree[idx][0]
            del self.cp_tree[child]
            del self.tree[idx]
        elif self.is_leaf(idx):
            parent = self.cp_tree[idx]
            del self.cp_tree[idx]
            if self.num_children(parent) == 1:
                del self.tree[parent]
            else:
                self.tree[parent].remove(idx)
        else:
            parent = self.cp_tree[idx]
            children = self.tree[idx]
            del self.cp_tree[idx]
            for child in children:
                self.cp_tree[child] = parent
                self.tree[parent].append(child)
                self.E[parent, child] = 1
            del self.tree[idx]

    def is_leaf(self, idx):
        if idx not in self.tree.keys():
            return True
        else:
            return False

    def is_root(self, idx):
        if idx in self.tree.keys() and idx not in self.cp_tree.keys():
            return True
        else:
            return False

    def num_children(self, idx):
        if self.is_leaf(idx):
            return 0
        else:
            return len(self.tree[idx])

# creates a readme file with the command in it. 
def write_readme(dname, args, script_name = os.path.basename(__file__)):
    readme_fname = dname + 'README.txt'
    open(readme_fname, 'w').close() # clear readme
    msg =  '    executed: ' + str(datetime.now()) + '\n'
    msg += 'command used:\n'
    msg += '\t```\n'
    msg += '\t' + ' '.join(['python', script_name] + [ '--' + str(k) + ' ' + _arg_val_to_str(v) for k, v in args.iteritems() ]) + '\n'
    msg += '\t```\n'
    fm.append_to_file(readme_fname, msg)
    readme = open(dname + "parameters.txt", 'w')
    for key, value in args.items():
        readme.write(str(key) + ":" + str(value) + "\n")
    readme.close()
    return readme_fname

def _arg_val_to_str(v):
    if isinstance(v, list):
        return ' '.join([ str(x) for x in v ])
    return str(v)

#  input: F (np.array) [m, l+r] mixed copy number of l breakpoints, r segments across m samples
#         Q (np.array) [l, r] binary indicator that breakpoint is in segment
#         num_seg_subsamples (int) number of segments (in addition to those containing breakpoints)
#             that are to be randomly kept in F
# output: F (np.array) [m, l+r'] r' is reduced number of segments
#         Q (np.array) [l, r']
#         org_indices (list of int) for each segment in output, the index of where it is found in input F
def randomly_remove_segments(F_phasing, Q, Q_unsampled, num_seg_subsamples):
    #print(Q)
    if num_seg_subsamples is None:
        return F_phasing, Q, Q_unsampled, None
    l_g, r = Q.shape
    l_g, r = int(l_g), int(r)
    g_un = Q_unsampled.shape[0]

    bp_segs = []
    for s in xrange(0, r):
        if sum(Q[:, s]): # segment s has a breakpoint in it
            bp_segs.append(s)
    for s in xrange(0, r):
        if sum(Q_unsampled[:, s]): # segment s has a breakpoint in it
            bp_segs.append(s)
    non_bp_segs = [ s for s in xrange(0, r) if s not in bp_segs ]  # all non breakpoint containing segments
    num_seg_subsamples = min(num_seg_subsamples, len(non_bp_segs)) # ensure not removing more segs than we have
    if num_seg_subsamples == len(non_bp_segs):
        return F_phasing, Q, Q_unsampled, None

    keeps = random_subset(non_bp_segs, num_seg_subsamples) # segments to keep
    keeps = set(sorted(bp_segs + keeps))
    drops = [ s for s in xrange(0, r) if s not in keeps ]

    Q = np.delete(Q, drops, axis = 1) # remove columns for segments we do not keep
    Q_unsampled = np.delete(Q_unsampled, drops, axis=1)
    #F = np.delete(F, [ s + l_g for s in drops ], axis = 1)
    F_phasing = np.delete(F_phasing, [ s + l_g + r for s in drops ], axis=1)
    F_phasing = np.delete(F_phasing, [s + l_g for s in drops], axis=1)
    # F_info_phasing = np.delete(F_info_phasing, [ s + l_g + r for s in drops ], axis=1)
    # F_info_phasing = np.delete(F_info_phasing, [s + l_g for s in drops], axis=1)

    return F_phasing, Q, Q_unsampled, [ s + l_g for s in keeps]

# returns a subset of lst containing k random elements
def random_subset(lst, k):
    result = []
    n = 0
    for item in lst:
        n += 1
        if len(result) < k:
            result.append(item)
        else:
            s = int(random.random() * n)
            if s < k:
                result[s] = item
    return result

def setup_get_UCE(args):
    return sv.get_UCE(*args)

def printnow(s):
    sys.stdout.write(s)
    sys.stdout.flush()


# # # # # # # # # # # # # # # #
#   W R I T E   O U T P U T   #
# # # # # # # # # # # # # # # #

#  input: F (np.array) [m, l+r] mixed copy number for all l bps and r segments for each sample
#         C (np.array) [n, l+r] integer copy number for each of n clones for all l bps and r' subset of r segments
#         org_indices (list of int) for each segment in F, the index of where it is found in input F_all
#         G (np.array) [l, l] G[i, j] == G[j, i] == 1 iff breakpoint i and j are mates. 0 otherwise
#         bp_attr (dict) key is breakpoint index. val is tuple (chrm (str), pos (int), extends_left (bool))
#         cv_attr (dict) key (int) is segment index. val is tuple (chrm (str), bgn_pos (int), end_pos (int))
# output: w (vcf_help.Writer) writer to be used to write entire .vcf file
def build_vcf_writer(F_phasing_full, C, org_indices, G, Q, bp_attr, cv_attr, metadata_fname):
    print(org_indices)
    m, l_g_2r = F_phasing_full.shape
    n, l_g_2rp = C.shape
    l, _ = G.shape
    g_2r = l_g_2r - l
    l_g = Q.shape[0]
    r = (l_g_2r - l_g)/2
    g = l_g - l
    print(C[:,:].shape, g_2r)

    if org_indices is not None: # only fill in values for segments not used if did not use some segments
        org_indices_minor = [org_indices[i] + r for i in range(len(org_indices))]
        c_org_indices = [ i for i in xrange(0, l_g) ] + org_indices + org_indices_minor
        print(c_org_indices, len(c_org_indices))
        C_out = -1*np.ones((n, l_g+2*r), dtype = float) # C with segments that were removed inserted back in with avg from F_full
        print(C_out.shape, C_out[:, c_org_indices].shape)
        C_out[:, c_org_indices] = C[:, :]           #   -1 is an indicator that this column should be omitted in validation
        C = C_out

    w = vh.Writer(m, n, metadata_fname)
    bp_ids = np.array([ 'bp' + str(b+1) for b in xrange(0, l) ], dtype = STR_DTYPE)
    for b in xrange(0, l): # force a breakpoint to not be mated with self
        G[b, b] = 0
    for b in xrange(0, l):
        chrm, pos, ext_left = bp_attr[b]
        rec_id = bp_ids[b]
        mate_id = bp_ids[np.where(G[b, :])[0][0]]
        fs = list(F_phasing_full[:, b])
        cps = list(C[:, b])
        if cps[0] < 0:
            cps = []
        w.add_bp(chrm, pos, ext_left, rec_id, mate_id, fs, cps)
    snv_ids = [ 'snv' + str(s+1) for s in xrange(0, g) ]

    cv_ids = [ 'cnv' + str(s+1) for s in xrange(0, r) ]
    for s in xrange(0, r):
        chrm, bgn, end = cv_attr[s]
        rec_id = cv_ids[s]
        fs = list(F_phasing_full[:, s + l_g])
        cps = list(C[:, s + l_g])
        if cps[0] < 0:
            cps = []
        w.add_cv(chrm, bgn, end, rec_id, fs, cps)

    return w

# d (str) is local directory path. all others are np.array
# input: F (np.array) [m, l+r'] mixed copy number for l bps, r' subset of r segments for each of m samples
#        F_full (np.array) [m, l+r] mixed copy number for all l bps and r segments for each sample
#        org_indices (list of int) for each segment in F, the index of where it is found in input F_all
#        writer (vcf_help.Writer) writer to be used to write entire .vcf file
def write_to_files(d, l_g, U, C, E, R, W, W_SV, W_SNV, W_SNV_UNSAMPLED, W_con, obj_val, F_phasing_full, F_unsampled_phasing_full, org_indices, writer, E_pre, R_pre, W_pre, B, A):
    l_g_2r = F_phasing_full.shape[1]
    r = (l_g_2r - l_g)/2
    n, _ = C.shape

    if org_indices is not None:
        org_indices_minor = [org_indices[i] + r for i in range(len(org_indices))]
        c_org_indices = [ i for i in xrange(0, l_g) ] + org_indices + org_indices_minor
        C_out = -1*np.ones((n, l_g+2*r), dtype = float) # C with segments that were removed inserted back in with avg from F_full
        C_out[:, c_org_indices] = C[:, :]           #   -1 is an indicator that this column should be omitted in validation
    else:
        C_out = C

    fnames = [ d + fname for fname in ['U.tsv', 'C.tsv', 'T.dot', 'F.tsv',  'W.tsv', 'obj_val.txt', 'unmixed.vcf', 'unmixed.xml','F_phasing_full.tsv','F_unsampled_phasing_full.tsv', 'W_SV.tsv', 'W_SNV_sampled.tsv', 'W_SNV_unsampled.tsv', 'W_CONCATENATE.tsv', 'T_pre.dot', 'B.tsv', 'A.tsv'] ]
    for fname in fnames:
        fm.touch(fname)
    np.savetxt(fnames[0], U, delimiter = '\t', fmt = '%.8f')
    np.savetxt(fnames[1], C_out, delimiter = '\t', fmt = '%.8f')
    np.savetxt(fnames[4], W, delimiter = '\t', fmt = '%d')
    np.savetxt(fnames[10], W_SV, delimiter='\t', fmt='%d')
    np.savetxt(fnames[11], W_SNV, delimiter='\t', fmt='%d')
    np.savetxt(fnames[12], W_SNV_UNSAMPLED, delimiter='\t', fmt='%d')
    np.savetxt(fnames[13], W_con, delimiter='\t', fmt='%d')
    np.savetxt(fnames[14], B, delimiter='\t', fmt='%d')
    np.savetxt(fnames[15], A, delimiter='\t', fmt='%d')
    np.savetxt(fnames[5], np.array([obj_val]), delimiter = '\t', fmt = '%.8f')
    np.savetxt(fnames[8], F_phasing_full, delimiter='\t', fmt='%.8f')
    np.savetxt(fnames[9], F_unsampled_phasing_full, delimiter='\t', fmt='%.8f')
    #writer.write(open(fnames[6], 'w'))
    dot = to_dot(E, R, W)
    open(fnames[2], 'w').write(dot.source) # write tree T in dot format
    dot.format = 'svg'
    dot.render(d + 'T')                    # display tree T in .svg
    dot = to_dot(E_pre, R_pre, W_pre)
    open(fnames[14], 'w').write(dot.source)  # write tree T in dot format
    dot.format = 'svg'
    dot.render(d + 'T_pre')
    write_xml(fnames[7], E, C, l_g)

#  input: E (np.array of int) [2n-1, 2n-1] 0 if no edge, 1 if edge between nodes i and j
#         R (np.array of int) [2n-1, 2n-1] cost of each edge in the tree
#         W (np.array of int) [2n-1, l] W[i, b] == 1 iff breakpoint b appears at node v_i. 0 otherwise
# output: dot (graphviz.dot.Digraph) directed tree representation of E
def to_dot(E, R, W):
    N = len(E)
    dot = Digraph(format = 'png')
    dot.node(str(N-1))
    for i in xrange(N-1, -1, -1):
        for j in xrange(N-1, -1, -1):
            if int(E[i, j]) == 1:
                num_breakpoints = sum(W[j, :])
                edge_label = ' ' + str(int(R[i, j])) + '/' + str(num_breakpoints)
                dot.node(str(j))
                dot.edge(str(i), str(j), label = edge_label)
    return dot

#  input: E (np.array)
def write_xml(fname, E, C, l_g):
    n, _ = E.shape

    root = Tree()
    root.name = str(n - 1)
    stack = [root]
    while stack:
        cur = stack.pop()
        i = int(cur.name)
        child_idxs = np.where(E[i, :] == 1)[0]
        for ci in child_idxs:
            child = cur.add_child(name = str(ci))
            child.dist = np.linalg.norm( np.subtract( C[i, l_g:], C[ci, l_g:] ), ord = 1 )
            stack.append(child)

    newick_str = root.write(features = ['name'], format = 1, format_root_node = True) # format_root_node=True puts root node name in str
    newick_tree = Phylo.read(StringIO(newick_str), 'newick') # format=1 gives branch lengths and names for all nodes (leaves and internal)

    for clade in newick_tree.find_clades():
        if clade.confidence is not None: # Phylo.read() stupidly interprets names of internal nodes as confidences for newick strings
            clade.name = clade.confidence
            clade.confidence = None
    xmltree = newick_tree.as_phyloxml() # convert to PhyloXML.Phylogeny type
    Phylo.write(xmltree, open(fname, 'w'), 'phyloxml')


# # # # # # # # # # # # # # # # # # # #
#   I N P U T   V A L I D A T I O N   #
# # # # # # # # # # # # # # # # # # # #

# input: Q (np.array of 0 or 1) [l, r] q_b,s == 1 if breakpoint b is in segment s. 0 otherwise
#        G (np.array of 0 or 1) [l, l] g_s,t == 1 if breakpoints s and t are mates. 0 otherwise
#        A (np.array of int) [m, l] a_p,b is number of mated reads for breakpoint b in sample p
#        H (np.array of int) [m, l] h_p,b is number of total reads for breakpoint b in sample p
#  does: exits with error message if any of the input is not valid
def check_valid_input(Q, Q_unsampled, G, A, H,F_phasing_full, F_unsampled_phasing_full):  ### A and H are empty matrices
    print("check valid input")
    l_g, r = np.shape(Q)
    print(l_g, r)
    l, _ = np.shape(G)
    g = l_g - l
    m = np.shape(A)[0]
    Q_msg = 'There is an issue with input binary matrix Q (indicates which segment each breakpoint belongs to). Each breakpoint must belong to exactly one segment.'
    Q_unsampled_msg = 'There is an issue with input binary matrix Q (indicates which segment each SNV belongs to). Each SNV must belong to exactly one segment.'

    G_msg = 'There is an issue with input binary matrix G (indicates which breakpoints are mates). Each breakpoint must be mated into pairs.'
    A_msg = 'There is an issue with input integer matricies A and H (indicating the number of reads mapped to each mated breakpoint and the number of total reads mapping to a breakpoint). The number of mated reads must be less or equal to the total reads and both should be non negative.'
    print(Q[np.where(np.sum(Q, axis=1) != 1)])
    sys.stdout.flush()

    raiseif(not np.all(np.sum(Q, 1) == 1), Q_msg)
    raiseif(not np.all(np.sum(Q_unsampled, 1) == 1), Q_unsampled_msg)

    print(np.where(np.sum(G, 0) != 2), np.where(np.sum(G, 0) != 2))
    raiseif(not np.all(np.sum(G, 0) == 2) or not np.all(np.sum(G, 1) == 2), G_msg)
    for i in xrange(0, l):
        for j in xrange(0, l):
            raiseif(G[i, j] != G[j, i], G_msg)
            raiseif(i == j and G[i, j] != 1, G_msg)

    for p in xrange(0, m):
        for b in xrange(0, l):
            raiseif(A[p, b] < 0 or A[p, b] > H[p, b], A_msg)
    return Q, Q_unsampled, G, A, H, F_phasing_full, F_unsampled_phasing_full

# raises exception if boolean is true
def raiseif(should_raise, msg):
    if should_raise:
        raise Exception(msg)

    # condition for G and A and H

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
    parser = argparse.ArgumentParser(prog = 'tusv.py', description = "unmixes mixed copy numbers for breakpoints and segments and infers phylogeny with various phylogenetic constraints")
    parser.add_argument('-i', '--input_directory', required = True, type = lambda x: fm.valid_dir_ext(parser, x, '.vcf'), help = 'directory containing a .vcf for each sample from a single patient')
    parser.add_argument('-o', '--output_directory', required = True, type = lambda x: fm.valid_dir(parser, x), help = 'empty directory for output U.tsv, C.tsv, and T.dot files to go')
    set_non_dir_args(parser)
    return vars(parser.parse_args(argv))

def set_non_dir_args(parser):
    parser.add_argument('-n', '--num_leaves', required = True, type = lambda x: fm.valid_int_in_range(parser, x, 2, MAX_NUM_LEAVES), help = 'number of leaves for inferred binary tree. total number of nodes will be 2*n-1')
    parser.add_argument('-c', '--c_max', required = True, type = lambda x: fm.valid_int_in_range(parser, x, 1, MAX_COPY_NUM), help = 'maximum allowed copy number at any node in the tree')
    parser.add_argument('-l', '--lambda1', default = 0.25, type = lambda x: fm.valid_float_above(parser, x, 0.0), help = 'regularization term to weight total tree cost against unmixing error in objective function. setting as 0.0 will put no tree cost constraint. setting as 1.0 will equally consider tree cost and unmixing error.')
    parser.add_argument('-a', '--lambda2', default = 6.25, type = lambda x: fm.valid_float_above(parser, x, 0.0), help = 'regularization term to weight error in inferred ratio between copy number of a breakpoint and the copy number of the segment originally containing the position of breakpoint')
    parser.add_argument('-t', '--cord_desc_iters', required = True, type = lambda x: fm.valid_int_in_range(parser, x, 1, MAX_CORD_DESC_ITERS), help = 'maximum number of cordinate descent iterations for each initialization of U')
    parser.add_argument('-r', '--restart_iters', required = True, type = lambda x: fm.valid_int_in_range(parser, x, 1, MAX_RESTART_ITERS), help = 'number of random initializations for picking usage matrix U')
    parser.add_argument('-p', '--processors', default = 1, type = lambda x: fm.valid_int_in_range(parser, x, 1, NUM_CORES), help = 'number of processors to use')
    parser.add_argument('-m', '--time_limit', type = int, help = 'maximum time (in seconds) allowed for a single iteration of the cordinate descent algorithm')
    parser.add_argument('-s', '--num_subsamples', type = int, default = None, help = 'number of segments (in addition to those containing breakpoints) that are to be randomly kept for deconvolution. default keeps all segments.')
    parser.add_argument('-d', '--metadata_file', default = METADATA_FNAME, type = lambda x: fm.is_valid_file(parser, x), help = 'file containing metadata information for output .vcf file')
    parser.add_argument('-b', '--overide_lambdas', action = 'store_true', help = 'specify this argument if you would like the parameters lambda1 and lambda2 to be set proportional to the input data set')
    parser.add_argument('-C', '--constant', default = 120, type = int, help = 'scaling constant for sampling SNVs')
    parser.add_argument('-sv_ub', '--sv_upperbound', default = -1, type = int, help = 'scaling constant for sampling SVs')
    parser.add_argument('-leaf', '--only_leaf', action = 'store_true', help = 'if only deconvolute for leaves')
    parser.add_argument('-col', '--collapse', action='store_true', help='if collapse nodes')
    parser.add_argument('-th', '--threshold', default = 0.0, type = lambda x: fm.valid_float_above(parser, x, 0.0), help = 'mean frequency threshold to collapsing')
    parser.add_argument('-scan', '--multi_num_clones', action='store_true', help='Scan a range of number of clones to get optimal number of clones')

# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == "__main__":
    main(sys.argv[1:])
