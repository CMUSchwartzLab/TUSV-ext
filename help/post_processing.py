#     file: file_manager.py
#   author: Jesse Eaton
#  created: 10/15/2017
# modified: 10/15/2017
#  purpose: Manages files used in arguments


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import numpy as np

from ete3 import Tree # for calculatin robinson foulds distance

from os.path import isfile, join


# # # # # # # # # #
#   P U B L I C   #
# # # # # # # # # #

#  input: dirname (str) name of directory where FNAMES are located
# output: C (np.array of int) [N, l+r] c_k,s is copy number of segment or breakpoint s in clone k
#         U (np.array of float) [m, N] u_p,k is percent of clone k making up sample p
#         T (ete3.coretype.tree.TreeNode) directed tree representing phylogeny
def get_CUT(dirname):
	C = np.genfromtxt(dirname + 'C.tsv', dtype = float)
	U = np.genfromtxt(dirname + 'U.tsv', dtype = float)
	T = _get_T(open(dirname + 'T.dot'))
	return C, U, T


# # # # # # # # # # #
#   P R I V A T E   #
# # # # # # # # # # #

def _get_T(file):
	txt = file.read()
	lines = txt.split('\n')
	lines.pop(0)
	lines.pop(len(lines)-1)
	lines = [ line.replace('\t', '') for line in lines ]
	lines = [ line.replace(' ', '') for line in lines ]

	edges = []
	edge_labels = {} # key is child node label. val is edge label
	for line in lines:
		if '->' in line:
			i, j = line.split('->')
			if '[' in j:
				j, edge_label = j.split('[')
				edges.append((i, j))
				edge_label = '[' + edge_label
				edge_labels[j] = edge_label
			else:
				edges.append((i, j))
				edge_labels[j] = ''

	parent_child_table = []
	for i, j in edges:
		if i in edge_labels.keys():
			i = i + edge_labels[i]
		if j in edge_labels.keys():
			j = j + edge_labels[j]
		parent_child_table.append((i, j, 1))
	return Tree.from_parent_child_table(parent_child_table)
