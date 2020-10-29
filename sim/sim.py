
#     file: sim.py
#   author: Jingyi Wang
#  created: 10/13/2017
#  modified: 10/13/2017
#  input: m, n, lambda, mu, outdir
#  output: u.tsv, c.tsv, t.dot, sample.vcf

##########
# Import #
##########
from __future__ import division
import sys
import os
import argparse 
import vcf
import itertools
import random
import operator
import datetime
import shutil

import numpy as np
import graphviz as gv
import chrm_prof as chpr
import gene_prof as gnpr

sys.path.insert(0, 'helper/')

import combine_copy_nums as ccn

#############
# Functions #
#############

def printnow(s, newline = True):
    s = str(s)
    if newline:
        s += '\n'
    sys.stdout.write(s)
    sys.stdout.flush()


def main(argv):
	args = get_args(argv)
	
	# input arguments
	m = args['m']
	n = args['n']
	num_mutes = args['num_mutes']
	
	directory = os.path.dirname(os.path.realpath(__file__))

	# for num_mutes in [10, 50, 100]:
		# for m in [1, 3, 5, 10]:
			# for n in [2, 3, 4, 5]:
				# print 'n:', n, 'm:', m, 'num_mutes:', num_mutes

	size_mutes = args['size_mutes']
	metaFile = args['meta_file']
	output_folder = args['output_folder']

	constants_dict = dict()
	constants_dict['mut_types'] = ['amp', 'rem', 'inv']
	constants_dict['exp_mut_size'] = size_mutes # default exp_mut_size is 5745000
	constants_dict['exp_mut_count'] = num_mutes / ( 2 * n - 2)
	constants_dict['cov'] = 20
	constants_dict['read_len'] = 300
	constants_dict['num_patients'] = 5
	
	# remove chrom_dict later
	chrom_dict = dict()
	chrom_dict[('1', 0)] = chpr.ChrmProf(248956422)
	chrom_dict[('1', 1)] = chpr.ChrmProf(248956422)
	chrom_dict[('2', 0)] = chpr.ChrmProf(242193529)
	chrom_dict[('2', 1)] = chpr.ChrmProf(242193529)
	chrom_dict[('3', 0)] = chpr.ChrmProf(198295559)
	chrom_dict[('3', 1)] = chpr.ChrmProf(198295559)

	# sub_folder_name = 'n_' + str(n) + '_m_' + str(m) + '_l_' + str(num_mutes)
	for patient_idx in range(1, 1 + constants_dict['num_patients']):
		patient_folder_name = 'patient' + str(patient_idx)

		# outputFolder = directory + '/sim_data' + '/' + sub_folder_name + '/' + patient_folder_name
		outputFolder = output_folder + '/' + patient_folder_name

		# clean up existing files under outputFolder
		if os.path.exists(outputFolder):
			shutil.rmtree(outputFolder)
		os.makedirs(outputFolder)
	
		l = random_get_tree(n) # list
		edge_list = get_edges(l)
		
		gp = gnpr.GeneProf(chrom_dict, constants_dict)

		t = Tree(edge_list, gp)

		geneprof_list = list()

		t.add_mutations_along_edges(t.rootNode, geneprof_list)

		generate_t(t, 'T.dot', outputFolder)

		U = random_get_usages(m, 2 * n - 1)

		l, sv_cn_idx_dict = get_bp_copy_num_idx_dict(t, n, constants_dict)
		r, seg_cn_idx_dict, seg_bgn_idx_dict, seg_end_idx_dict = get_seg_copy_num_idx_dict(t, n)
		C = generate_c(t, n, constants_dict)

		c_p, c_m = generate_seg_cp_paternal(t, n)

		F = generate_f(U, C)

		a, h, mate_dict = get_a_h_mate_dict(t, n, constants_dict)

		output_tsv(U, '/U.tsv', outputFolder)
		output_tsv(C, '/C.tsv', outputFolder)
		output_tsv(F, '/F.tsv', outputFolder)

		generate_s(metaFile, t, l, sv_cn_idx_dict, r, seg_cn_idx_dict, seg_bgn_idx_dict, seg_end_idx_dict, F, U, C, c_p, c_m, a, h, mate_dict, outputFolder)


# given a number n, generate all possible directed binary trees with n nodes.
# eg. if n = 4, return [ [1,[1,[1,1]]], [1,[[1,1],1]], [[1,1],[1,1]], [[1,[1,1]],1], [[[1,1],1],1] ]
# each 1 in the list represents a tree node.
def all_possible_trees(n):
	l = list()
	for tree in all_possible_trees_helper(n):
		l.append(tree)
	return l


def all_possible_trees_helper(n):
	# base case
	if n == 1:
		yield 1

	# recursive case
	for i in range(1, n):
		left_list = all_possible_trees_helper(i)
		right_list = all_possible_trees_helper(n - i)
		for left, right in itertools.product(left_list, right_list):
			yield [left, right]


# input n (number of leaves), output a random tree represented by a list.
# tree node is represented by 1.
def random_get_tree(n):
	l = all_possible_trees(n)
	idx = random.randint(0, len(l) - 1)
	return l[idx]


# input a tree (list), return number of leaves
def get_number_of_leaves(l):
	if l[0] == 1 and l[1] == 1:
		return 2
	elif l[0] == 1 and l[1] != 1:
		return 1 + get_number_of_leaves(l[1])
	elif [0] != 1 and l[1] == 1:
		return 1 + get_number_of_leaves(l[0])
	else:
		return get_number_of_leaves(l[0]) + get_number_of_leaves(l[1])


# given a tree (list), return list of edges(tuple). tuple format: (id, parent_id, l/r)
# l/r: left or right child
def get_edges(l):
	n = get_number_of_leaves(l)
	p = 2 * n - 1 # idx of root node
	leaf_list = list(range(1, n + 1)) # indices of leaf nodes
	result = list()
	get_edges_helper(l, n, leaf_list, p, result)
	return result


def get_edges_helper(l, n, leaf_list, p, result):
	left, right = l[0], l[1]
	if left == 1 and right == 1:
		result.append((leaf_list[0], p, 'l'))
		result.append((leaf_list[1], p, 'r'))
	elif left == 1 and right != 1:
		result.append((leaf_list[0], p, 'l'))
		result.append((p - 1, p, 'r'))
		get_edges_helper(l[1], n, leaf_list[1:], p - 1, result)
	elif right == 1 and left != 1:
		result.append((leaf_list[-1], p, 'r'))
		result.append((p-1, p, 'l'))
		get_edges_helper(l[0], n-1, leaf_list[:-1], p - 1, result)
	else:
		n_left = get_number_of_leaves(l[0])
		leaf_list_left = list(range(1, n_left + 1))
		leaf_list_right = list(range(n_left + 1, n + 1))
		result.append((p - 2, p, 'l')) # left
		result.append((p - 1, p, 'r')) # right
		get_edges_helper(l[0], n_left, leaf_list_left, p - 2, result) # left
		get_edges_helper(l[1], n, leaf_list_right, p - 1, result) # right
	return result


# pwd: os.path.dirname(os.path.realpath(__file__))
# generate a dot file to save the random tree with n nodes and random mutations
def generate_t(tree, filename, directory):
	dot = gv.Digraph()
	nodes = tree.node_list
	edges = get_graph_edges(tree.edge_list)
	add_nodes(dot, nodes)
	add_edges(dot, edges)
	# print dot.source
	dot.save(filename, directory)
	return 


def add_nodes(graph, nodes):
    for n in nodes:
        if isinstance(n, tuple):
            graph.node(str(n[0]), **n[1])
        else:
            graph.node(str(n))
    return graph


def add_edges(graph, edges):
    for e in edges:
        if isinstance(e[0], tuple):
            graph.edge(*e[0], **e[1])
        else:
            graph.edge(*e)
    return graph


# edge_list: list of (node_id (int), parent_node_id (int)) tuples. 
# return list of (from (str), to (str)) tuples.
def get_graph_edges(edge_list):
	result = list()
	for(node_id, parent_node_id, lr) in edge_list:
		result.append((str(parent_node_id), str(node_id)))
	return result

# return a 1 by n np array sum to one
def get_usage_for_one_patient(n):
	a = np.random.dirichlet(np.ones(n),size = 1)
	return a


# return m by n np array, each row sum to one
def random_get_usages(m, n):
	a = get_usage_for_one_patient(n)
	for i in range(1, m):
		b = get_usage_for_one_patient(n)
		a = np.concatenate((a, b), axis = 0)
	return a


# input a tree (Tree object) and n (number of leaf nodes)
# output a l (number of bps) and a dictionary
# key: chrom
# val: dictionary
#      key: (pos, isLeft) tuple. When two bps have the same pos, the one with isLeft == True comes first.
#      val: index of the bp
def get_bp_copy_num_idx_dict(tree, n, constants_dict):
	# put all bps in leaf nodes into dictionary d1
	d1 = dict()
	for idx in range(1, n + 1):
		temp_bp_dict = tree.idx_node_dict[idx].geneProf.get_sv_read_nums_dict(constants_dict['cov'], constants_dict['read_len'])
		for chrom in temp_bp_dict.keys():
			if chrom not in d1:
				d1[chrom] = set()
			for (pos, isLeft) in temp_bp_dict[chrom]:
				if (pos, isLeft) not in d1[chrom]:
					d1[chrom].add((pos, isLeft))

	# sort vals in d1 based on pos and isLeft (isLeft == True comes first) and chrom
	idx = 0
	sorted_chrom = sorted(d1.keys())
	d2 = dict()
	for chrom in sorted_chrom:
		d2[chrom] = dict()
		sorted_pos_list = sorted(list(set(map(operator.itemgetter(0), d1[chrom]))))
		for pos in sorted_pos_list:
			if (pos, True) in d1[chrom]:
				d2[chrom][(pos, True)] = idx
				idx += 1
			if (pos, False) in d1[chrom]:
				d2[chrom][(pos, False)] = idx
				idx += 1

	return idx, d2


# input a tree (Tree object) and n (number of leaf nodes)
# output a r and a dictionary
# key: chrom
# val: dictionary
#      key: (bgn, end) tuple
#      val: index
def get_seg_copy_num_idx_dict(tree, n):
	d1 = dict()
	for idx in range(1, n + 1):
		temp_copy_nums_dict = tree.idx_node_dict[idx].geneProf.get_copy_nums_dict()
		for chrom in temp_copy_nums_dict.keys():
			if chrom not in d1:
				d1[chrom] = list()
			(bgns, ends, cps) = temp_copy_nums_dict[chrom]
			d1[chrom].append([bgns, ends, cps])

	idx = 0
	sorted_chrom = sorted(d1.keys())
	d2 = dict()
	bgn2idx = dict()
	end2idx = dict()
	for chrom in sorted_chrom:
		d2[chrom] = dict()
		# random generate usages to match the input format
		usages = [float(1/len(d1[chrom]))] * len(d1[chrom])
		[res_bgns, res_ends, res_cps] = ccn.combine_copy_nums(d1[chrom], usages)

		for i in range(len(res_bgns)):
			d2[chrom][(res_bgns[i], res_ends[i])] = idx
			bgn2idx[(chrom, res_bgns[i])] = idx
			end2idx[(chrom, res_ends[i])] = idx
			idx += 1

	return idx, d2, bgn2idx, end2idx


def make_2d_list(rows, cols):
	result = list()
	for r in range(rows):
		result.append([0] * cols)
	return result


# loop through each node in tree(Tree), 
# for each treeNode: use self.get_copy_nums_dict() to get bgns, ends, cps list for each chromosomes
#                    use self.get_sv_read_nums_dict(cov, read_len) to get bps and their corresponding information for each chromosomes
# output c ((2n-1)*(l+r) matrix) 
def generate_c(tree, n, constants_dict):

	l, sv_cn_idx_dict = get_bp_copy_num_idx_dict(tree, n, constants_dict)
	r, seg_cn_idx_dict, seg_bgn_idx_dict, seg_end_idx_dict = get_seg_copy_num_idx_dict(tree, n)

	c = make_2d_list(len(tree.node_list), (l + r))
	for idx in tree.node_list:
		row = idx - 1

		# add copy number for break points
		temp_bp_dict = tree.idx_node_dict[idx].geneProf.get_sv_read_nums_dict(constants_dict['cov'], constants_dict['read_len'])
		for chrom in temp_bp_dict:
			for (pos, isLeft) in temp_bp_dict[chrom]:
				cp = temp_bp_dict[chrom][(pos, isLeft)]["copy_num"]
				col = sv_cn_idx_dict[chrom][(pos, isLeft)]
				c[row][col] = cp

		# add copy number for segments
		temp_copy_nums_dict = tree.idx_node_dict[idx].geneProf.get_copy_nums_dict()
		for chrom in temp_copy_nums_dict:
			(bgns, ends, cps) = temp_copy_nums_dict[chrom]
			for i in range(len(bgns)):
				cp = cps[i]
				seg_indices_list = get_indices_for_segment(seg_bgn_idx_dict, seg_end_idx_dict, (chrom, bgns[i]), (chrom, ends[i]))
				for j in range(len(seg_indices_list)):
					col = seg_indices_list[j] + l
					c[row][col] = cp

	result = np.array(c)

	return result


# given start and end position, output list of segment indices (continuous)
# s = (chrom, bgn_pos), e = (chrom, end_pos)
def get_indices_for_segment(bgn2idx, end2idx, s, e):
	result = list()
	firstIdx = bgn2idx[s]
	lastIdx = end2idx[e]
	for i in range(firstIdx, lastIdx + 1):
		result.append(i)
	return result


# return a ((2n-1) * r) matrix contains paternal chrom copy number for each segment
# and a ((2n-1) * r) matrix contains maternal chrom copy number for each segment
def generate_seg_cp_paternal(tree, n):
	r, seg_cn_idx_dict, seg_bgn_idx_dict, seg_end_idx_dict = get_seg_copy_num_idx_dict(tree, n)
	c_p = make_2d_list(len(tree.node_list), r)
	c_m = make_2d_list(len(tree.node_list), r)
	for idx in tree.node_list:
		row = idx - 1
		temp_chrom_dict = tree.idx_node_dict[idx].geneProf.chrom_dict
		# temp_copy_nums_dict = tree.idx_node_dict[idx].geneProf.get_copy_nums_dict()
		for chrom in list(filter(lambda x: x[1] == 0, temp_chrom_dict.keys())):
			(bgns, ends, cps) = temp_chrom_dict[chrom].get_copy_nums()
			for i in range(len(bgns)):
				cp = cps[i]
				seg_indices_list = get_indices_for_segment(seg_bgn_idx_dict, seg_end_idx_dict, (chrom[0], bgns[i]), (chrom[0], ends[i]))
				for col in seg_indices_list:
					c_p[row][col] = cp

		for chrom in list(filter(lambda x: x[1] == 1, temp_chrom_dict.keys())):
			(bgns, ends, cps) = temp_chrom_dict[chrom].get_copy_nums()
			for i in range(len(bgns)):
				cp = cps[i]
				seg_indices_list = get_indices_for_segment(seg_bgn_idx_dict, seg_end_idx_dict, (chrom[0], bgns[i]), (chrom[0], ends[i]))
				for col in seg_indices_list:
					c_m[row][col] = cp
	return c_p, c_m


# given u (m * (2n-1) matrix) and c ((2n-1)*(l+r) matrix), output f (m * (l+r) matrix)
def generate_f(u, c):
	return np.dot(u, c)


# return matrix a, matrix h, and dictionary mate_dict
# matrix a: a ((2n-1) * l) matrix contains mated_reads info
# matrix h: a ((2n-1) * l) matrix contains total_reads info
# mate_dict:
#           key: (chrom, pos, isLeft)
#           val: (mate_chrom, mate_pos, mate_isLeft)
def get_a_h_mate_dict(tree, n, constants_dict):
	l, sv_cn_idx_dict = get_bp_copy_num_idx_dict(tree, n, constants_dict)
	a, h = make_2d_list(len(tree.node_list), l), make_2d_list(len(tree.node_list), l)

	mate_dict = {}
	for node_name in tree.node_list:
		gene_prof = tree.idx_node_dict[node_name].geneProf
		sv_dict = gene_prof.get_sv_read_nums_dict(constants_dict['cov'], constants_dict['read_len'])
		for chrm in sv_dict.keys():
			for cur_pos, cur_is_left in sv_dict[chrm]:
				mat_pos, mat_is_left = sv_dict[chrm][(cur_pos, cur_is_left)]['mate']
				cur_key, mat_key = (chrm, cur_pos, cur_is_left), (chrm, mat_pos, mat_is_left)
				if cur_key not in mate_dict:
					mate_dict[cur_key] = mat_key
				elif mate_dict[cur_key] != mat_key:
					print 'There was an error generating SVs. Rerun sim.py until there are no errors.'
					print 'cur_key:\t' + str(cur_key)
					print 'mat_key:\t' + str(mat_key)
					print 'cur_key was already mated with:\t' + str(mate_dict[cur_key])
					exit()
				if mat_key not in mate_dict:
					mate_dict[mat_key] = cur_key
				elif mate_dict[mat_key] != cur_key:
					print 'There was an error generating SVs. Rerun sim.py until there are no errors.'
					print 'cur_key:\t' + str(cur_key)
					print 'mat_key:\t' + str(mat_key)
					print 'mat_key was already mated with:\t' + str(mate_dict[mat_key])
					exit()

				j = sv_cn_idx_dict[chrm][(cur_pos, cur_is_left)]
				a[node_name - 1][j] = sv_dict[chrm][(cur_pos, cur_is_left)]['mated_reads']
				h[node_name - 1][j] = sv_dict[chrm][(cur_pos, cur_is_left)]['total_reads']
	return a, h, mate_dict

# given a matrix, save as tsv file
def output_tsv(mtx, output_file, output_folder):
	with open(output_folder + output_file, "w") as f:
		f.write("\n".join("\t".join(map(str, x)) for x in mtx))


# given cnv_idx(int) and r(int), output rec_id(str)
# eg. given rec_idx = 1, r = 12, output 'cnv01'
def get_cnv_rec_id(cnv_idx, r):
	return 'cnv' + str(cnv_idx).zfill(len(str(r)))


# given sv_idx(int) and r(int), output rec_id(str)
def get_sv_rec_id(sv_idx, l):
	return 'sv' + str(sv_idx).zfill(len(str(l)))


# a, h, mate_dict = get_a_h_mate_dict(t, n, constants_dict)
# generate a vcf file for each sample
def generate_s(metaFile, tree, l, sv_cn_idx_dict, r, seg_cn_idx_dict, seg_bgn_idx_dict, seg_end_idx_dict, F, U, C, c_p, c_m, a, h, mate_dict, outputFolder):
	vcf_reader = vcf.Reader(open(metaFile, 'r'))
	vcf_reader.metadata['filedate'][0] = datetime.datetime.now().date().strftime('%Y%m%d') # set date to current date
	f_p = np.dot(U, c_p)
	f_m = np.dot(U, c_m)
	mixed_a = np.dot(U, a) # m * l
	mixed_h = np.dot(U, h) # m * l
	for i in range(len(U)):
		sample_idx = i + 1
		temp_file = outputFolder + '/sample' + str(sample_idx) + '.vcf'
		temp_writer = vcf.Writer(open(temp_file, 'w'), vcf_reader)
		alt_type, gt_cnv = 'CNV', '1|1' # constants for all cnv records
		for chrom in sorted(seg_cn_idx_dict.keys()):
			for (key, val) in sorted(seg_cn_idx_dict[chrom].items(), key = lambda x: x[1]):
				pos = key[0]
				rec_id = get_cnv_rec_id(val, r)
				info_end = key[1]
				cn = [f_p[i][val], f_m[i][val]]
				temp_writer.write_record(generate_cnv(chrom, pos, rec_id, alt_type, info_end, gt_cnv, cn))

		alt_ori, alt_cS, alt_wMA, gt_sv = True, str(), True, '1|0' # constants for all sv records
		for chrom in sorted(sv_cn_idx_dict.keys()):
			for (key, val) in sorted(sv_cn_idx_dict[chrom].items(), key = lambda x: x[1]):
				pos, isLeft = key[0], key[1]
				rec_id = get_sv_rec_id(val, l)
				(mate_chrom, mate_pos, mate_isLeft) = mate_dict[(chrom, pos, isLeft)]
				mate_id = sv_cn_idx_dict[mate_chrom][(mate_pos, mate_isLeft)]
				alt_chr, alt_pos = mate_chrom, mate_pos
				cnadj = F[i][val]
				bdp, dp = int(round(mixed_a[i][val])), int(round(mixed_h[i][val]))
				info_mateid = get_sv_rec_id(mate_id, l)
				alt_rO = False if mate_isLeft == True else True
				temp_writer.write_record(generate_sv(chrom, pos, rec_id, alt_chr, alt_pos, alt_ori, alt_rO, alt_cS, alt_wMA, info_mateid, gt_sv, cnadj, bdp, dp))


# chrom(str), pos(int), rec_id(str), ref(str), qual = None, filter(list), fmt = 'GT:CNADJ', sample = ['TUMOR', 'NORMAL']
# type(alts): list
#              type(alts[0]) = class 'vcf.model._Breakend'
#              dir(alts[0]): [..., 'chr', 'connectingSequence', 'orientation', 'pos', 'remoteOrientation', 'type', 'withinMainAssembly']
#                           eg. alts[0] = ]1:149965077]
#                               'chr' = str(1), 'connectingSequence' = str(), 'orientation' = True, 'pos' = int(149965077)
#                               'remoteOrientation' = False, 'type' = 'BND', 'withinMainAssembly' = True
# info(dict), info['SVTYPE'] = 'BND', info['MATEID'] = mate_sv_rec_id (str)
def generate_sv(chrm, pos, rec_id, alt_chr, alt_pos, alt_ori, alt_rO, alt_cS, alt_wMA, info_mateid, gt, cnadj, bdp, dp):
	ref = '.'
	alts = list()
	alts.append(vcf.parser._Breakend(alt_chr, alt_pos, alt_ori, alt_rO, alt_cS, alt_wMA))
	qual = None
	filt = list()
	info = dict()
	info['SVTYPE'] = 'BND'
	info['MATEID'] = info_mateid
	fmt = 'GT:CNADJ:BDP:DP'
	samples = ['TUMOR', 'NORMAL']
	calls = [vcf.model._Call(0, 'TUMOR', svCallData(gt,cnadj,bdp,dp)), vcf.model._Call(1, 'NORMAL', svCallData('0|0',0, 0, dp))]
	newRec = vcf.model._Record(chrm, pos, rec_id, ref, alts, qual, filt, info, fmt, samples, calls)
	return newRec


def generate_cnv(chrm, pos, rec_id, alt_type, info_end, gt, cn):
	ref = '.'
	alts = list()
	alts.append(vcf.model._SV(alt_type))
	qual = None
	filt = list()
	info = dict()
	info['IMPRECISE'] = True
	info['END'] = info_end
	fmt = 'GT:CN'
	samples = ['TUMOR', 'NORMAL']
	calls = [vcf.model._Call(0, 'TUMOR', cnvCallData(gt, cn)), vcf.model._Call(1, 'NORMAL', cnvCallData('0|0',[1, 1]))]
	newRec = vcf.model._Record(chrm, pos, rec_id, ref, alts, qual, filt, info, fmt, samples, calls)
	return newRec


def is_cnv_record(rec):
	return rec.ID[0:3] == 'cnv'


def is_sv_record(rec):
	return rec.ID[0:2] == 'sv'


#########
# Class #
#########

class Tree:

	def __init__(self, edge_list, gp):
		self.edge_list = edge_list
		self.geneProf = gp # without any mutation
		self.node_list = self.get_node_list()
		self.rootNode, self.idx_node_dict = self.construct_tree()


	def get_node_list(self):
		l = list()
		for (idx, parent_idx, lr) in self.edge_list:
			if idx not in l:
				l.append(idx)
			if parent_idx not in l:
				l.append(parent_idx)
		return l


	# construct the tree structure, return root node (treeNode) and a dictionary
	# key: treeNode.index, val: treeNode
	def construct_tree(self):
		d = dict()
		for (idx, parent_idx, lr) in self.edge_list:
			# add new node to d
			if parent_idx not in d:
				parentNode = TreeNode(parent_idx, self.geneProf)
				d[parent_idx] = parentNode
			if idx not in d:
				childNode = TreeNode(idx, self.geneProf)
				d[idx] = childNode
			# add left or right child
			if lr == 'l': # left child
				d[parent_idx].left = d[idx]
			elif lr == 'r': # right child
				d[parent_idx].right = d[idx]
			# add parent node
			d[idx].parent = d[parent_idx]
		for nodeIdx in d:
			if d[nodeIdx].parent == None: # root node
				rootNode = d[nodeIdx]
				return rootNode, d


	def print_tree_info(self):
		for idx in self.idx_node_dict:
			print 'node:', idx
			print self.idx_node_dict[idx].geneProf.print_chrm_seq()


	# print current node index, parent node index, left child index, and right child index
	def print_node_relation(self):
		for idx in self.idx_node_dict:
			print 'node:', self.idx_node_dict[idx].index
			if self.idx_node_dict[idx].parent != None:
				print 'parent node:', self.idx_node_dict[idx].parent.index 
			else:
				print 'no parent node!'
			if self.idx_node_dict[idx].left != None:
				print 'left node:', self.idx_node_dict[idx].left.index
			else:
				print 'no left child node!'
			if self.idx_node_dict[idx].right != None:
				print 'right node:', self.idx_node_dict[idx].right.index
			else:
				print 'no right child node!'
			print ""


	def get_number_of_leaves(self, l):
		if l[0] == 1 and l[1] == 1:
			return 2
		elif l[0] == 1 and l[1] != 1:
			return 1 + get_number_of_leaves(l[1])
		elif [0] != 1 and l[1] == 1:
			return 1 + get_number_of_leaves(l[0])
		else:
			return get_number_of_leaves(l[0]) + get_number_of_leaves(l[1])


	def print_node_info(self):
		for idx in self.idx_node_dict:
			print 'node', idx, ':', self.idx_node_dict[idx].geneProf.print_info()


	def print_node_gp(self):
		for idx in self.idx_node_dict:
			print 'node', idx, ':', self.idx_node_dict[idx].geneProf.print_chrm_seq()


	def add_mutations_along_edges(self, node, geneprof_list):
		if not node:
			return
		curr_gp = node.geneProf
		geneprof_list.append(curr_gp)
		# print 'node:', node.index, 'geneprof_list:', geneprof_list

		if node.left != None:
			curr_gp_copied_left = curr_gp.deepcopy()
			# reset copied_node.geneProf.mutCount and copied_node.geneProf.maxCount
			curr_gp_copied_left.mutCount, curr_gp_copied_left.maxCount = 0, curr_gp_copied_left.get_mut_count()

			curr_gp_copied_left.multi_mutations(geneprof_list) 
			node.left.geneProf = curr_gp_copied_left
			self.add_mutations_along_edges(node.left, geneprof_list)

		if node.right != None:
			curr_gp_copied_right = curr_gp.deepcopy()

			# reset copied_node.geneProf.mutCount and copied_node.geneProf.maxCount
			curr_gp_copied_right.mutCount, curr_gp_copied_right.maxCount = 0, curr_gp_copied_right.get_mut_count()

			curr_gp_copied_right.multi_mutations(geneprof_list)
			node.right.geneProf = curr_gp_copied_right
			self.add_mutations_along_edges(node.right, geneprof_list)
		return


class TreeNode:
	def __init__(self, index, gp):
		self.index = index
		self.geneProf = gp
		self.left = None
		self.right = None
		self.parent = None


class svCallData:
	def __init__(self, gt = '0|0', cnadj = '0', bdp = '0', dp = '100'):
		self.GT = gt
		self.CNADJ = str(cnadj)
		self.BDP = str(bdp)
		self.DP = str(dp)
		self.__getitem__ = self

	def __call__(self, var):
		return [self.CNADJ, self.BDP, self.DP]


class cnvCallData:
	def __init__(self, gt = '0|0', cns = [1, 1]):
		self.GT = gt
		self.CN = cns
		self.__getitem__ = self

	def __call__(self, var):
		return [','.join(str(x) for x in self.CN)]


###########################################
##### COMMAND LINE ARGUMENT FUNCTIONS #####
###########################################

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'sim.py', description = "generate U.tsv, C.tsv, T.dot, sample.vcf")
	parser.add_argument('-f', '--metadata_file', type = str, dest = "meta_file", required = True)
	parser.add_argument('-m', '--num_samples', type = int, dest = "m", required = True)
	parser.add_argument('-n', '--num_leaves', type = int, dest = "n", required = True)
	parser.add_argument('-c', '--total_number_of_mutations', type = int, dest = "num_mutes", required = True)
	parser.add_argument('-s', '--expect_mut_len', type = int, dest = "size_mutes", required = True)
	parser.add_argument('-o', '--output_folder', type = str, dest = "output_folder", required = True)
	return vars(parser.parse_args(argv))


##############################
##### CALL MAIN FUNCTION #####
##############################

if __name__ == "__main__":
	main(sys.argv[1:])