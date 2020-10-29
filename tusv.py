#     file: tusv.py
#   author: Jesse Eaton
#  created: 10/13/2017
# modified: 10/14/2017
#  purpose: Unmixes mixed copy numbers for breakpoints and segments and infers phylogeny
#             with various phylogenetic constraints


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

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

# custom modules
sys.path.insert(0, 'model/')
sys.path.insert(0, 'help/')
import solver as sv
import file_manager as fm      # sanitizes file and directory arguments
import generate_matrices as gm # gets F, Q, G, A, H from .vcf files
import printer as pt
import vcf_help as vh


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
	unmix(args['input_directory'], args['output_directory'], args['num_leaves'], args['c_max'], args['lambda1'], args['lambda2'], args['restart_iters'], args['cord_desc_iters'], args['processors'], args['time_limit'], args['metadata_file'], args['num_subsamples'], args['overide_lambdas'])

#  input: num_seg_subsamples (int or None) number of segments to include in deconvolution. these are
#           in addition to any segments contining an SV as thos are manditory for the SV. None is all segments
def unmix(in_dir, out_dir, n, c_max, lamb1, lamb2, num_restarts, num_cd_iters, num_processors, time_limit, metadata_fname, num_seg_subsamples, should_overide_lambdas):

	F_full, Q, G, A, H, bp_attr, cv_attr = gm.get_mats(in_dir)
	check_valid_input(Q, G, A, H)

	F, Q, org_indxs = randomly_remove_segments(F_full, Q, num_seg_subsamples)

	# replace lambda1 and lambda2 with input derived values if should_orveride_lamdas was specified
	if should_overide_lambdas:
		m = len(F)
		l, r = Q.shape
		lamb1 = float(l + r) / float(r) * float(m) / float(2 * (n-1) )
		lamb2 = float(l + r) / float(l)

	Us, Cs, Es, obj_vals, Rs, Ws = [], [], [], [], [], []
	num_complete = 0
	for i in xrange(0, num_restarts):
		U, C, E, R, W, obj_val, err_msg = sv.get_UCE(F, Q, G, A, H, n, c_max, lamb1, lamb2, num_cd_iters, time_limit)
		printnow(str(i + 1) + ' of ' + str(num_restarts) + ' random restarts complete\n')
		Us.append(U)
		Cs.append(C)
		Es.append(E)
		Rs.append(R)
		Ws.append(W)
		obj_vals.append(obj_val)

	best_i = 0
	best_obj_val = obj_vals[best_i]
	for i, obj_val in enumerate(obj_vals):
		if obj_val < best_obj_val:
			best_obj_val = obj_val
			best_i = i

	writer = build_vcf_writer(F_full, Cs[best_i], org_indxs, G, bp_attr, cv_attr, metadata_fname)

	write_to_files(out_dir, Us[best_i], Cs[best_i], Es[best_i], Rs[best_i], Ws[best_i], F, obj_vals[best_i], F_full, org_indxs, writer)

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
def randomly_remove_segments(F, Q, num_seg_subsamples):
	if num_seg_subsamples is None:
		return F, Q, None
	l, r = Q.shape
	l, r = int(l), int(r)

	bp_segs = []
	for s in xrange(0, r):
		if sum(Q[:, s]): # segment s has a breakpoint in it
			bp_segs.append(s)
	non_bp_segs = [ s for s in xrange(0, r) if s not in bp_segs ]  # all non breakpoint containing segments
	num_seg_subsamples = min(num_seg_subsamples, len(non_bp_segs)) # ensure not removing more segs than we have
	if num_seg_subsamples == len(non_bp_segs):
		return F, Q, None

	keeps = random_subset(non_bp_segs, num_seg_subsamples) # segments to keep
	keeps = sorted(bp_segs + keeps)
	drops = [ s for s in xrange(0, r) if s not in keeps ]

	Q = np.delete(Q, drops, axis = 1) # remove columns for segments we do not keep
	F = np.delete(F, [ s + l for s in drops ], axis = 1)
	
	return F, Q, [ s + l for s in keeps ]

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
def build_vcf_writer(F, C, org_indices, G, bp_attr, cv_attr, metadata_fname):
	m, _ = F.shape
	n, _ = C.shape
	l, _ = G.shape
	r = F.shape[1] - l
	
	if org_indices is not None: # only fill in values for segments not used if did not use some segments
		c_org_indices = [ i for i in xrange(0, l) ] + org_indices
		C_out = -1*np.ones((n, l+r), dtype = float) # C with segments that were removed inserted back in with avg from F_full
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
		fs = list(F[:, b])
		cps = list(C[:, b])
		if cps[0] < 0:
			cps = []
		w.add_bp(chrm, pos, ext_left, rec_id, mate_id, fs, cps)
	cv_ids = [ 'cnv' + str(s+1) for s in xrange(0, r) ]
	for s in xrange(0, r):
		chrm, bgn, end = cv_attr[s]
		rec_id = cv_ids[s]
		fs = list(F[:, s + l])
		cps = list(C[:, s + l])
		if cps[0] < 0:
			cps = []
		w.add_cv(chrm, bgn, end, rec_id, fs, cps)

	return w

# d (str) is local directory path. all others are np.array
# input: F (np.array) [m, l+r'] mixed copy number for l bps, r' subset of r segments for each of m samples
#        F_full (np.array) [m, l+r] mixed copy number for all l bps and r segments for each sample
#        org_indices (list of int) for each segment in F, the index of where it is found in input F_all
#        writer (vcf_help.Writer) writer to be used to write entire .vcf file
def write_to_files(d, U, C, E, R, W, F, obj_val, F_full, org_indices, writer):
	l = F.shape[1]
	if org_indices is not None:
		l = F.shape[1] - len(org_indices)
	r = F_full.shape[1] - l
	n, _ = C.shape
	
	if org_indices is not None:
		c_org_indices = [ i for i in xrange(0, l) ] + org_indices
		C_out = -1*np.ones((n, l+r), dtype = float) # C with segments that were removed inserted back in with avg from F_full
		C_out[:, c_org_indices] = C[:, :]           #   -1 is an indicator that this column should be omitted in validation
	else:
		C_out = C

	fnames = [ d + fname for fname in ['U.tsv', 'C.tsv', 'T.dot', 'F.tsv', 'W.tsv', 'obj_val.txt', 'unmixed.vcf', 'unmixed.xml'] ]
	for fname in fnames:
		fm.touch(fname)
	np.savetxt(fnames[0], U, delimiter = '\t', fmt = '%.8f')
	np.savetxt(fnames[1], C_out, delimiter = '\t', fmt = '%.8f')
	np.savetxt(fnames[3], F_full, delimiter = '\t', fmt = '%.8f')
	np.savetxt(fnames[4], W, delimiter = '\t', fmt = '%.8f')
	np.savetxt(fnames[5], np.array([obj_val]), delimiter = '\t', fmt = '%.8f')
	writer.write(open(fnames[6], 'w'))
	dot = to_dot(E, R, W)
	open(fnames[2], 'w').write(dot.source) # write tree T in dot format
	dot.format = 'svg'
	dot.render(d + 'T')                    # display tree T in .svg
	write_xml(fnames[7], E, C, l)

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
def write_xml(fname, E, C, l):
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
			child.dist = np.linalg.norm( np.subtract( C[i, l:], C[ci, l:] ), ord = 1 )
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
def check_valid_input(Q, G, A, H):
	l, r = np.shape(Q)
	m = np.shape(A)[0]
	Q_msg = 'There is an issue with input binary matrix Q (indicates which segment each breakpoint belongs to). Each breakpoint must belong to exactly one segment.'
	G_msg = 'There is an issue with input binary matrix G (indicates which breakpoints are mates). Each breakpoint must be mated into pairs.'
	A_msg = 'There is an issue with input integer matricies A and H (indicating the number of reads mapped to each mated breakpoint and the number of total reads mapping to a breakpoint). The number of mated reads must be less or equal to the total reads and both should be non negative.'

	raiseif(not np.all(np.sum(Q, 1) == 1), Q_msg)

	raiseif(not np.all(np.sum(G, 0) == 2) or not np.all(np.sum(G, 0) == 2), G_msg)
	for i in xrange(0, l):
		for j in xrange(0, l):
			raiseif(G[i, j] != G[j, i], G_msg)
			raiseif(i == j and G[i, i] != 1, G_msg)

	for p in xrange(0, m):
		for b in xrange(0, l):
			raiseif(A[p, b] < 0 or A[p, b] > H[p, b], A_msg)

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

# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
	main(sys.argv[1:])
