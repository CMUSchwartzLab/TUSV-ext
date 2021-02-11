
#     file: sim.py
#   author: Jingyi Wang
#  created: 10/13/2017
#  modified: 10/13/2017
#  input: m, n, lambda, mu, outdir
#  output: u.tsv, c.tsv, t.dot, sample.vcf

#  co-author: Xuecong Fu
#  modified: 2021.1-3

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
import pickle

import numpy as np
import graphviz as gv
import chrm_prof as chpr
import gene_prof as gnpr
from chrm_prof import *


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
	num_mutes_snv = args['num_mutes_snv']
	num_patients = args['num_patients']
	
	directory = os.path.dirname(os.path.realpath(__file__))

	# for num_mutes in [10, 50, 100]:
		# for m in [1, 3, 5, 10]:
			# for n in [2, 3, 4, 5]:
				# print 'n:', n, 'm:', m, 'num_mutes:', num_mutes

	size_mutes = args['size_mutes']
	metaFile = args['meta_file']
	output_folder = args['output_folder']

	constants_dict = dict()
	constants_dict['mut_types'] = ['amp', 'inv', 'rem','trans']
	constants_dict['exp_mut_size'] = size_mutes # default exp_mut_size is 5745000
	constants_dict['exp_mut_count'] = num_mutes / ( 2 * n - 2)
	if num_mutes_snv is not None:
		constants_dict['snv_mut_lambda'] = num_mutes_snv / (2*n - 2)
	else:
		constants_dict['snv_mut_lambda'] = None
	constants_dict['cov'] = 20
	constants_dict['read_len'] = 300
	constants_dict['num_patients'] = num_patients
	constants_dict['num_leaves'] = n
	constants_dict['num_samples'] = m
	
	# remove chrom_dict later
	chrom_dict = dict()
	chrom_dict[('1', 0)] = chpr.ChrmProf(1000, '1', 0)
	chrom_dict[('1', 1)] = chpr.ChrmProf(1000, '1', 1)
	# chrom_dict[('2', 0)] = chpr.ChrmProf(1000, '2', 0)
	# chrom_dict[('2', 1)] = chpr.ChrmProf(1000, '2', 1)
	# chrom_dict[('3', 0)] = chpr.ChrmProf(198295559, '3', 0)
	# chrom_dict[('3', 1)] = chpr.ChrmProf(198295559, '3', 1)

	# chrom_dict[('1', 0)] = chpr.ChrmProf(249198692, '1', 0) #from ICGC CNV data
	# chrom_dict[('1', 1)] = chpr.ChrmProf(249198692, '1', 1)
	# chrom_dict[('2', 0)] = chpr.ChrmProf(243048760, '2', 0)
	# chrom_dict[('2', 1)] = chpr.ChrmProf(243048760, '2', 1)
	# chrom_dict[('3', 0)] = chpr.ChrmProf(197856433, '3', 0)
	# chrom_dict[('3', 1)] = chpr.ChrmProf(197856433, '3', 1)
	# chrom_dict[('4', 0)] = chpr.ChrmProf(190921709, '4', 0)
	# chrom_dict[('4', 1)] = chpr.ChrmProf(190921709, '4', 1)
	# chrom_dict[('5', 0)] = chpr.ChrmProf(186092833, '5', 0)
	# chrom_dict[('5', 1)] = chpr.ChrmProf(186092833, '5', 1)
	# chrom_dict[('6', 0)] = chpr.ChrmProf(170918031, '6', 0)
	# chrom_dict[('6', 1)] = chpr.ChrmProf(170918031, '6', 1)
	# chrom_dict[('7', 0)] = chpr.ChrmProf(159119220, '7', 0)
	# chrom_dict[('7', 1)] = chpr.ChrmProf(159119220, '7', 1)
	# chrom_dict[('8', 0)] = chpr.ChrmProf(146293414, '8', 0)
	# chrom_dict[('8', 1)] = chpr.ChrmProf(146293414, '8', 1)
	# chrom_dict[('9', 0)] = chpr.ChrmProf(141071475, '9', 0)
	# chrom_dict[('9', 1)] = chpr.ChrmProf(141071475, '9', 1)
	# chrom_dict[('10', 0)] = chpr.ChrmProf(135434551, '10', 0)
	# chrom_dict[('10', 1)] = chpr.ChrmProf(135434551, '10', 1)
	# chrom_dict[('11', 0)] = chpr.ChrmProf(134944770, '11', 0)
	# chrom_dict[('11', 1)] = chpr.ChrmProf(134944770, '11', 1)
	# chrom_dict[('12', 0)] = chpr.ChrmProf(133777645, '12', 0)
	# chrom_dict[('12', 1)] = chpr.ChrmProf(133777645, '12', 1)
	# chrom_dict[('13', 0)] = chpr.ChrmProf(115106996, '13', 0)
	# chrom_dict[('13', 1)] = chpr.ChrmProf(115106996, '13', 1)
	# chrom_dict[('14', 0)] = chpr.ChrmProf(107285437, '14', 0)
	# chrom_dict[('14', 1)] = chpr.ChrmProf(107285437, '14', 1)
	# chrom_dict[('15', 0)] = chpr.ChrmProf(102400037, '15', 0)
	# chrom_dict[('15', 1)] = chpr.ChrmProf(102400037, '15', 1)
	# chrom_dict[('16', 0)] = chpr.ChrmProf(90163275, '16', 0)
	# chrom_dict[('16', 1)] = chpr.ChrmProf(90163275, '16', 1)
	# chrom_dict[('17', 0)] = chpr.ChrmProf(81048659, '17', 0)
	# chrom_dict[('17', 1)] = chpr.ChrmProf(81048659, '17', 1)
	# chrom_dict[('18', 0)] = chpr.ChrmProf(78015057, '18', 0)
	# chrom_dict[('18', 1)] = chpr.ChrmProf(78015057, '18', 1)
	# chrom_dict[('19', 0)] = chpr.ChrmProf(59095126, '19', 0)
	# chrom_dict[('19', 1)] = chpr.ChrmProf(59095126, '19', 1)
	# chrom_dict[('20', 0)] = chpr.ChrmProf(62912463, '20', 0)
	# chrom_dict[('20', 1)] = chpr.ChrmProf(62912463, '20', 1)
	# chrom_dict[('21', 0)] = chpr.ChrmProf(48084820, '21', 0)
	# chrom_dict[('21', 1)] = chpr.ChrmProf(48084820, '21', 1)
	# chrom_dict[('22', 0)] = chpr.ChrmProf(51219006, '22', 0)
	# chrom_dict[('22', 1)] = chpr.ChrmProf(51219006, '22', 1)
	# chrom_dict[('23', 0)] = chpr.ChrmProf(155233846, '23', 0)
	# chrom_dict[('23', 1)] = chpr.ChrmProf(155233846, '23', 1)

	# sub_folder_name = 'n_' + str(n) + '_m_' + str(m) + '_l_' + str(num_mutes)
	if not os.path.exists(output_folder):
		os.mkdir(output_folder)
	readme = open(output_folder + "/README.md", 'w')
	for key, value in constants_dict.items():
		readme.write(str(key) + ":" + str(value) + "\n")
	readme.close()
	for patient_idx in range(1, 1 + constants_dict['num_patients']):
		patient_folder_name = 'patient' + str(patient_idx)

		# outputFolder = directory + '/sim_data' + '/' + sub_folder_name + '/' + patient_folder_name
		outputFolder = output_folder + '/' + patient_folder_name

		# clean up existing files under outputFolder
		if os.path.exists(outputFolder):
			shutil.rmtree(outputFolder)
		os.makedirs(outputFolder)
	
		l = random_get_tree(n) # list
		print(l)
		edge_list = get_edges(l)  ###xf: generate edges list with format of [(0,1,'r'/'l'),...]
		
		gp = gnpr.GeneProf(chrom_dict, constants_dict)

		t = Tree(edge_list, gp)
		print(t.node_list)

		geneprof_list = list()

		t.add_mutations_along_edges(t.rootNode, geneprof_list)

		generate_t(t, 'T.dot', outputFolder)

		U = random_get_usages(m, 2 * n - 1)

		l, sv_cn_idx_dict = get_bp_copy_num_idx_dict(t, n, constants_dict)
		r, seg_cn_idx_dict, seg_bgn_idx_dict, seg_end_idx_dict = get_seg_copy_num_idx_dict(t, n)
		### xf: combine the segment settings from both two alleles and also different node from mutations
		bool_list = np.random.choice([True, False], r)
		if constants_dict['snv_mut_lambda'] is None:
			C = generate_c(t, n, constants_dict, bool_list)
			c_p, c_m = generate_seg_cp_paternal(t, n, bool_list)
			F = generate_f(U, C)
			a, h, mate_dict = get_a_h_mate_dict(t, n, constants_dict)
			generate_s(metaFile, t, l, sv_cn_idx_dict, r, seg_cn_idx_dict, seg_bgn_idx_dict, seg_end_idx_dict, F, U, C,
					   c_p, c_m, a, h, mate_dict, outputFolder)
			output_tsv(U, '/U.tsv', outputFolder)
			output_tsv(C, '/C.tsv', outputFolder)
			output_tsv(F, '/F.tsv', outputFolder)

		else:
			g, snv_cn_idx_dict = get_snv_copy_num_idx_dict(t)
			C, C_unsampled_snv, snv_sampled_idx, snv_unsampled_idx = generate_c_snv(t, n, constants_dict, bool_list)
			c_p, c_m = generate_seg_cp_paternal(t, n, bool_list)
			F = generate_f(U, C)
			F_unsampled_snv = generate_f(U, C_unsampled_snv)
			a, h, mate_dict = get_a_h_mate_dict(t, n, constants_dict)
			generate_s_snv(metaFile, t, l, sv_cn_idx_dict, r, seg_cn_idx_dict, g, snv_cn_idx_dict, snv_sampled_idx,
					   snv_unsampled_idx, seg_bgn_idx_dict, seg_end_idx_dict, F,
					   F_unsampled_snv, U, C, c_p, c_m, a, h, mate_dict, outputFolder)
			output_tsv(U, '/U.tsv', outputFolder)
			output_tsv(C, '/C.tsv', outputFolder)
			output_tsv(F, '/F.tsv', outputFolder)
			output_tsv(F_unsampled_snv, '/F_unsampled_snv.tsv', outputFolder)
			output_tsv(C_unsampled_snv, '/C_unsampled_snv.tsv', outputFolder)

		edge_list_pickle = open(outputFolder + "/edge_list.pickle", 'wb')
		pickle.dump(edge_list, edge_list_pickle)
		edge_list_pickle.close()
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
def random_get_tree(n):  ###xf: select a random tree from all possible trees
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
	steiner_list = list(range(n, p, 1)) ### xf: debug the tree generation part by adding a fixed steiner node list
	result = list()
	get_edges_helper(l, n, leaf_list, p, result, steiner_list)
	return result


def get_edges_helper(l, n, leaf_list, p, result, steiner_node_list): ### xf: debug the tree generation part by adding a fixed steiner node list
	left, right = l[0], l[1]
	if left == 1 and right == 1:
		result.append((leaf_list[0], p, 'l'))
		result.append((leaf_list[1], p, 'r'))
	elif left == 1 and right != 1:
		result.append((leaf_list[0], p, 'l'))
		steiner_node = steiner_node_list.pop()
		result.append((steiner_node, p, 'r'))
		get_edges_helper(l[1], n-1, leaf_list[1:], steiner_node, result, steiner_node_list)
	elif right == 1 and left != 1:
		result.append((leaf_list[-1], p, 'r'))
		steiner_node = steiner_node_list.pop()
		result.append((steiner_node, p, 'l'))
		get_edges_helper(l[0], n-1, leaf_list[:-1], steiner_node, result, steiner_node_list)
	else:
		n_left = get_number_of_leaves(l[0])
		leaf_list_left = leaf_list[:n_left]
		leaf_list_right = leaf_list[n_left:]
		steiner_node1 = steiner_node_list.pop()
		steiner_node2 = steiner_node_list.pop()
		result.append((steiner_node1, p, 'l')) # left
		result.append((steiner_node2, p, 'r')) # right
		get_edges_helper(l[0], n_left, leaf_list_left, steiner_node1, result, steiner_node_list) # left
		get_edges_helper(l[1], n - n_left, leaf_list_right, steiner_node2, result, steiner_node_list) # right
	return result, steiner_node_list


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

### xf: add snv copy number idx dict
def get_snv_copy_num_idx_dict(tree):
	d = set()
	d2 = {}
	for idx in tree.node_list:
		temp_snv_dict = tree.idx_node_dict[idx].geneProf.get_snv_dict()

		for tuple in temp_snv_dict.keys():
				if tuple not in d:
					d.add(tuple)
	sorted_d = sorted(list(d),key=lambda x: (x[0], x[1]))
	for k in range(len(sorted_d)):
		d2[sorted_d[k]] = k

	return len(list(d)), d2



# input a tree (Tree object) and n (number of leaf nodes)
# output a l (number of bps) and a dictionary
# key: chrom
# val: dictionary
#      key: (pos, isLeft) tuple. When two bps have the same pos, the one with isLeft == True comes first.
#      val: index of the bp
def get_bp_copy_num_idx_dict(tree, n, constants_dict):
	# put all bps in leaf nodes into dictionary d1
	d1 = dict()
	d_chrom = {}
	for idx in tree.node_list:
		temp_bp_dict = tree.idx_node_dict[idx].geneProf.get_sv_read_nums_dict(constants_dict['cov'], constants_dict['read_len'])

		###xf: for all chromosomes for both alleles
		for chrom in temp_bp_dict.keys():
			if chrom not in d1:
				d1[chrom] = set()
			for (pos, isLeft, chr_) in temp_bp_dict[chrom]:
				if (pos, isLeft, chr_) not in d1[chrom]:
					d1[chrom].add((pos, isLeft, chr_))


	# sort vals in d1 based on pos and isLeft (isLeft == True comes first) and chrom
	idx = 0
	sorted_chrom = sorted(d1.keys())
	d2 = dict()
	for chrom in sorted_chrom:
		#print(d1[chrom])
		d2[chrom] = dict()
		# for i in range(0,len(d1[chrom])):
		# 	for j in range(i+1,len(d1[chrom])):
		# 		if list(d1[chrom])[i][0] == list(d1[chrom])[j][0]:
		# 			print(list(d1[chrom])[i])
		sorted_pos_list = sorted(list(set(map(operator.itemgetter(0), d1[chrom]))))
		#print(len(sorted_pos_list))
		for pos in sorted_pos_list:
			if (pos, True, chrom) in d1[chrom]:
				d2[chrom][(pos, True, chrom)] = idx
				idx += 1
			if (pos, False, chrom) in d1[chrom]:
				d2[chrom][(pos, False, chrom)] = idx
				idx += 1
		#print("d1", chrom, len(d1[chrom]))
		#print("d2",chrom, len(d2[chrom]))
	return idx, d2


# input a tree (Tree object) and n (number of leaf nodes)
# output a r and a dictionary
# key: chrom
# val: dictionary
#      key: (bgn, end) tuple
#      val: index
def get_seg_copy_num_idx_dict(tree, n):
	d1 = dict()
	for idx in tree.node_list: ### xf: go through all the nodes in the tree
		temp_copy_nums_dict = tree.idx_node_dict[idx].geneProf.get_copy_nums_dict()
		for chrom in temp_copy_nums_dict.keys():
			if chrom not in d1:
				d1[chrom] = list()
			(bgns, ends, cps1, cps2) = temp_copy_nums_dict[chrom]
			d1[chrom].append([bgns, ends, cps1, cps2])

	idx = 0
	sorted_chrom = sorted(d1.keys())
	d2 = dict()
	bgn2idx = dict()
	end2idx = dict()
	for chrom in sorted_chrom:
		d2[chrom] = dict()
		# random generate usages to match the input format
		usages = [float(1/len(d1[chrom]))] * len(d1[chrom])
		[res_bgns, res_ends, res_cps1, res_cps2] = ccn.combine_copy_nums_quartet(d1[chrom], usages)

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


def generate_c_snv(tree, n, constants_dict, bool_list, subsample=0.001):

	l, sv_cn_idx_dict = get_bp_copy_num_idx_dict(tree, n, constants_dict)
	r, seg_cn_idx_dict, seg_bgn_idx_dict, seg_end_idx_dict = get_seg_copy_num_idx_dict(tree, n)
	g, snv_cn_idx_dict = get_snv_copy_num_idx_dict(tree)

	g_sample = int(g*subsample)
	snv_sampled_idx = np.sort(np.random.choice(g, size=g_sample, replace=False))
	snv_unsampled_idx = np.setdiff1d(np.arange(g), snv_sampled_idx)

	c = make_2d_list(len(tree.node_list), (l + g_sample + 2*r))
	c_unsampled_snv = make_2d_list(len(tree.node_list), (g-g_sample))
	for idx in tree.node_list:
		row = idx - 1
		# add copy number for break points
		temp_bp_dict = tree.idx_node_dict[idx].geneProf.get_sv_read_nums_dict(constants_dict['cov'], constants_dict['read_len'])
		for chrom in temp_bp_dict:
			for (pos, isLeft, chr_,) in temp_bp_dict[chrom]:
				cp = temp_bp_dict[chrom][(pos, isLeft, chr_)]["copy_num"]
				col = sv_cn_idx_dict[chrom][(pos, isLeft, chr_)]
				c[row][col] = cp

		temp_snv_dict = tree.idx_node_dict[idx].geneProf.get_snv_dict()
		for (chrm, pos) in temp_snv_dict.keys():
			if snv_cn_idx_dict[(chrm, pos)] in snv_sampled_idx:
				cp = temp_snv_dict[(chrm, pos)]["copy_num"]
				col = np.where(snv_sampled_idx == snv_cn_idx_dict[(chrm, pos)])[0][0] + l
				c[row][col] = cp
			else:
				cp = temp_snv_dict[(chrm, pos)]["copy_num"]
				col = np.where(snv_unsampled_idx == snv_cn_idx_dict[(chrm, pos)])[0][0]
				c_unsampled_snv[row][col] = cp

		# add copy number for segments
		temp_copy_nums_dict = tree.idx_node_dict[idx].geneProf.get_copy_nums_dict()
		for chrom in temp_copy_nums_dict:
			(bgns, ends, cps1, cps2) = temp_copy_nums_dict[chrom]
			for i in range(len(bgns)):
				cp1 = cps1[i]
				cp2 = cps2[i]
				seg_indices_list = get_indices_for_segment(seg_bgn_idx_dict, seg_end_idx_dict, (chrom, bgns[i]), (chrom, ends[i]))
				for j in range(len(seg_indices_list)):
					col = seg_indices_list[j] + l + g_sample
					if bool_list[col-l-g_sample]:
						c[row][col] = cp1
						c[row][col + r] = cp2
					else:
						c[row][col] = cp2
						c[row][col + r] = cp1
	result = np.array(c)
	return result, np.array(c_unsampled_snv), snv_sampled_idx, snv_unsampled_idx

# loop through each node in tree(Tree), 
# for each treeNode: use self.get_copy_nums_dict() to get bgns, ends, cps list for each chromosomes
#                    use self.get_sv_read_nums_dict(cov, read_len) to get bps and their corresponding information for each chromosomes
# output c ((2n-1)*(l+r) matrix) ### xf: --> (2n-1)*(l+2r)
def generate_c(tree, n, constants_dict, bool_list):

	l, sv_cn_idx_dict = get_bp_copy_num_idx_dict(tree, n, constants_dict)
	r, seg_cn_idx_dict, seg_bgn_idx_dict, seg_end_idx_dict = get_seg_copy_num_idx_dict(tree, n)


	c = make_2d_list(len(tree.node_list), (l + 2*r))
	for idx in tree.node_list:
		row = idx - 1
		# add copy number for break points
		temp_bp_dict = tree.idx_node_dict[idx].geneProf.get_sv_read_nums_dict(constants_dict['cov'], constants_dict['read_len'])
		for chrom in temp_bp_dict:
			for (pos, isLeft, chr_,) in temp_bp_dict[chrom]:
				cp = temp_bp_dict[chrom][(pos, isLeft, chr_)]["copy_num"]
				col = sv_cn_idx_dict[chrom][(pos, isLeft, chr_)]
				c[row][col] = cp

		# add copy number for segments
		temp_copy_nums_dict = tree.idx_node_dict[idx].geneProf.get_copy_nums_dict()
		for chrom in temp_copy_nums_dict:
			(bgns, ends, cps1, cps2) = temp_copy_nums_dict[chrom]
			for i in range(len(bgns)):
				cp1 = cps1[i]
				cp2 = cps2[i]
				seg_indices_list = get_indices_for_segment(seg_bgn_idx_dict, seg_end_idx_dict, (chrom, bgns[i]),
														   (chrom, ends[i]))
				for j in range(len(seg_indices_list)):
					col = seg_indices_list[j] + l
					if bool_list[col - l]:
						c[row][col] = cp1
						c[row][col + r] = cp2
					else:
						c[row][col] = cp2
						c[row][col + r] = cp1

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
def generate_seg_cp_paternal(tree, n, bool_list):
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
					if bool_list[col]:
						c_p[row][col] = cp
					else:
						c_m[row][col] = cp

		for chrom in list(filter(lambda x: x[1] == 1, temp_chrom_dict.keys())):
			(bgns, ends, cps) = temp_chrom_dict[chrom].get_copy_nums()
			for i in range(len(bgns)):
				cp = cps[i]
				seg_indices_list = get_indices_for_segment(seg_bgn_idx_dict, seg_end_idx_dict, (chrom[0], bgns[i]), (chrom[0], ends[i]))
				for col in seg_indices_list:
					if bool_list[col]:
						c_m[row][col] = cp
					else:
						c_p[row][col] = cp
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
		#print("node",node_name)
		gene_prof = tree.idx_node_dict[node_name].geneProf
		sv_dict = gene_prof.get_sv_read_nums_dict(constants_dict['cov'], constants_dict['read_len'])
		for chrm in sv_dict.keys():
			#print("ch",chrm)
			for cur_pos, cur_is_left, cur_chr in sv_dict[chrm]:
				mat_pos, mat_is_left, mat_chr = sv_dict[chrm][(cur_pos, cur_is_left, cur_chr)]['mate']
				cur_key, mat_key = (cur_chr, cur_pos, cur_is_left), (mat_chr, mat_pos, mat_is_left)
				#print(cur_key, mat_key)
				if cur_key not in mate_dict:
					mate_dict[cur_key] = mat_key
				elif mate_dict[cur_key] != mat_key:
					#print(sv_dict[cur_chr])
					print 'There was an error generating SVs. Rerun sim.py until there are no errors.'
					print 'cur_key:\t' + str(cur_key)
					print 'mat_key:\t' + str(mat_key)
					print 'cur_key was already mated with:\t' + str(mate_dict[cur_key])
					exit()
				if mat_key not in mate_dict:
					mate_dict[mat_key] = cur_key
				elif mate_dict[mat_key] != cur_key:
					#print(sv_dict[mat_chr])
					print 'There was an error generating SVs. Rerun sim.py until there are no errors.'
					print 'cur_key:\t' + str(cur_key)
					print 'mat_key:\t' + str(mat_key)
					print 'mat_key was already mated with:\t' + str(mate_dict[mat_key])
					exit()

				j = sv_cn_idx_dict[chrm][(cur_pos, cur_is_left, cur_chr)]
				a[node_name - 1][j] = sv_dict[chrm][(cur_pos, cur_is_left, cur_chr)]['mated_reads']
				h[node_name - 1][j] = sv_dict[chrm][(cur_pos, cur_is_left, cur_chr)]['total_reads']
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

### xf: add for snv
def get_snv_rec_id(snv_idx, g):
	return 'snv' + str(snv_idx).zfill(len(str(g)))


# a, h, mate_dict = get_a_h_mate_dict(t, n, constants_dict)
# generate a vcf file for each sample
def generate_s_snv(metaFile, tree, l, sv_cn_idx_dict, r, seg_cn_idx_dict, g, snv_cn_idx_dict, snv_sampled_idx, snv_unsampled_idx,
			   seg_bgn_idx_dict, seg_end_idx_dict, F, F_unsampled_snv, U, C, c_p, c_m, a, h, mate_dict, outputFolder):
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
		temp_file_unsampled_snv = outputFolder + '/sample' + str(sample_idx) + '_unsampled_snv.vcf'
		temp_writer_unsampled_snv = vcf.Writer(open(temp_file_unsampled_snv, 'w'), vcf_reader)
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
				#print(sv_cn_idx_dict[mate_chrom])
				mate_id = sv_cn_idx_dict[mate_chrom][(mate_pos, mate_isLeft, mate_chrom)]
				alt_chr, alt_pos = mate_chrom, mate_pos
				cnadj = F[i][val]
				bdp, dp = int(round(mixed_a[i][val])), int(round(mixed_h[i][val]))
				info_mateid = get_sv_rec_id(mate_id, l)
				alt_rO = False if mate_isLeft == True else True
				temp_writer.write_record(generate_sv(chrom, pos, rec_id, alt_chr, alt_pos, alt_ori, alt_rO, alt_cS, alt_wMA, info_mateid, gt_sv, cnadj, bdp, dp))
		###xf: add snvs
		for (key, val) in sorted(snv_cn_idx_dict.items(), key= lambda x: x[1]):
			chrm, pos = key
			rec_id = get_snv_rec_id(val, g)
			gt_snv = '0|1'
			if val in snv_sampled_idx:
				cnadj_snv = F[i][np.where(snv_sampled_idx == val)[0][0]+l]
				temp_writer.write_record(generate_snv(chrm, pos, rec_id, gt_snv, cnadj_snv))
			else:
				cnadj_snv = F_unsampled_snv[i][np.where(snv_unsampled_idx == val)[0][0]]
				temp_writer_unsampled_snv.write_record(generate_snv(chrm, pos, rec_id, gt_snv, cnadj_snv))


def generate_s(metaFile, tree, l, sv_cn_idx_dict, r, seg_cn_idx_dict,
			   seg_bgn_idx_dict, seg_end_idx_dict, F, U, C, c_p, c_m, a, h, mate_dict, outputFolder):
	vcf_reader = vcf.Reader(open(metaFile, 'r'))
	vcf_reader.metadata['filedate'][0] = datetime.datetime.now().date().strftime('%Y%m%d')  # set date to current date
	f_p = np.dot(U, c_p)
	f_m = np.dot(U, c_m)
	mixed_a = np.dot(U, a)  # m * l
	mixed_h = np.dot(U, h)  # m * l
	for i in range(len(U)):
		sample_idx = i + 1
		temp_file = outputFolder + '/sample' + str(sample_idx) + '.vcf'
		temp_writer = vcf.Writer(open(temp_file, 'w'), vcf_reader)

		alt_type, gt_cnv = 'CNV', '1|1'  # constants for all cnv records
		for chrom in sorted(seg_cn_idx_dict.keys()):
			for (key, val) in sorted(seg_cn_idx_dict[chrom].items(), key=lambda x: x[1]):
				pos = key[0]
				rec_id = get_cnv_rec_id(val, r)
				info_end = key[1]
				cn = [f_p[i][val], f_m[i][val]]
				temp_writer.write_record(generate_cnv(chrom, pos, rec_id, alt_type, info_end, gt_cnv, cn))

		alt_ori, alt_cS, alt_wMA, gt_sv = True, str(), True, '1|0'  # constants for all sv records
		for chrom in sorted(sv_cn_idx_dict.keys()):
			for (key, val) in sorted(sv_cn_idx_dict[chrom].items(), key=lambda x: x[1]):
				pos, isLeft = key[0], key[1]
				rec_id = get_sv_rec_id(val, l)
				(mate_chrom, mate_pos, mate_isLeft) = mate_dict[(chrom, pos, isLeft)]
				# print(sv_cn_idx_dict[mate_chrom])
				mate_id = sv_cn_idx_dict[mate_chrom][(mate_pos, mate_isLeft, mate_chrom)]
				alt_chr, alt_pos = mate_chrom, mate_pos
				cnadj = F[i][val]
				bdp, dp = int(round(mixed_a[i][val])), int(round(mixed_h[i][val]))
				info_mateid = get_sv_rec_id(mate_id, l)
				alt_rO = False if mate_isLeft == True else True
				temp_writer.write_record(
					generate_sv(chrom, pos, rec_id, alt_chr, alt_pos, alt_ori, alt_rO, alt_cS, alt_wMA, info_mateid,
								gt_sv, cnadj, bdp, dp))



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

def generate_snv(chrm, pos, rec_id, gt, cnadj):
	ref = '.'
	alt = list()
	qual = None
	filt = list()
	info = dict()
	fmt = 'GT:CNADJ'
	samples = ['TUMOR', 'NORMAL']
	calls = [vcf.model._Call(0, 'TUMOR', snvCallData(gt,cnadj)), vcf.model._Call(1, 'NORMAL', snvCallData('0|0', 0))]
	newRec = vcf.model._Record(chrm, pos, rec_id, ref, alt, qual, filt, info, fmt, samples, calls)
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


def is_snv_record(rec):
	return rec.ID[0:3] == 'snv'


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


	def add_mutations_along_edges(self, node, geneprof_list): ### xf: node is Treenode class
		if not node:
			return
		curr_gp = node.geneProf
		geneprof_list.append(curr_gp) ### xf: make sure each calling of add_mutations_along_edges will be saved, node is resursive
		print 'node:', node.index, 'geneprof_list:'
		for k in curr_gp.chrom_dict.keys():
			print(k)
			c = curr_gp.chrom_dict[k].mut
			while c != None:
				print(c.bgn, c.end, chpr._get_org_pos(c, True)[0], chpr._get_org_pos(c, False)[0])
				c = c.r

		if node.left != None:
			curr_gp_copied_left = curr_gp.deepcopy()
			# reset copied_node.geneProf.mutCount and copied_node.geneProf.maxCount
			curr_gp_copied_left.mutCount, curr_gp_copied_left.maxCount = 0, curr_gp_copied_left.get_mut_count() ### xf: get_mut_count: random.poisson(exp_mut_rate)

			curr_gp_copied_left.multi_mutations(geneprof_list) ### xf: generate multiple mutations


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

class snvCallData:
	def __init__(self, gt = '0|0', cnadj = '0'):
		self.GT = gt
		self.CNADJ = str(cnadj)
		self.__getitem__ = self

	def __call__(self, var):
		return [self.CNADJ]


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
	parser.add_argument('-cs', '--total_number_of_mutations_snv', type=int, dest="num_mutes_snv", required=True)
	parser.add_argument('-s', '--expect_mut_len', type = int, dest = "size_mutes", required = True)
	parser.add_argument('-o', '--output_folder', type = str, dest = "output_folder", required = True)
	parser.add_argument('-p', '--num_patients', type=int, dest="num_patients", default=5)
	return vars(parser.parse_args(argv))


##############################
##### CALL MAIN FUNCTION #####
##############################

if __name__ == "__main__":
	main(sys.argv[1:])