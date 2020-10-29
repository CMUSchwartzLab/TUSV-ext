#     file: analyze_results.py
#   author: Jesse Eaton
#  created: 1/28/2018
# modified: 1/28/2018
#  purpose: analyzing results from TCGA data run through tusv for purpose of differentiating
#             recurring versus non recurring tumors


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import re       # for extracting string from between quotations
import numpy as np
import scipy as sp

import matplotlib.pyplot as plt
import matplotlib
font = {'family' : 'normal',
        'size'   : 20}
matplotlib.rc('font', **font)

sys.path.insert(0, '../help/')
import file_manager as fm   # sanitizes file and directory arguments
import post_processing as pp


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

FNAMES = ['F.tsv', 'U.tsv', 'C.tsv', 'W.tsv', 'T.dot', 'obj_val.txt', 'unmixed.vcf', 'unmixed.xml']


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)

	rec_dir, nrc_dir = args['recur_directory'], args['non_recur_directory']
	# rec_subdir_names, nrc_subdir_names = sorted(fm.get_subdir_names(rec_dir)), sorted(fm.get_subdir_names(nrc_dir))
	nrc_subdir_names = sorted(fm.get_subdir_names(nrc_dir))

	# rec_stats = _get_stats(rec_dir, rec_subdir_names)
	nrc_stats = _get_stats(nrc_dir, nrc_subdir_names)

	for key in nrc_stats:

		# rec_stat = rec_stats[key]
		nrc_stat = nrc_stats[key]

		printnow(key)
		# for group_name, stat in [ ('non-recurring', nrc_stat), ('recurring', rec_stat) ]:
		for group_name, stat in [ ('non-recurring', nrc_stat) ]:
			printnow('\t' + group_name + ':')
			printnow('\t\tmean:\t' + str(np.mean(stat)))
			printnow('\t\tstdev:\t' + str(np.std(stat)))
		printnow('')

	colors = ['#2ecc71', '#3498db', '#9b59b6']
	color = colors[2]

	for key in nrc_stats:
		parts = plt.violinplot(nrc_stats[key])
		# plt.violinplot(rec_stats[key], [2])
		if key == 'num_unique_clones':
			plt.ylabel('Number of Unique Clones')
		elif key == 'num_bps':
			plt.ylabel('Number of Breakpoints')
		elif key == 'diversities':
			plt.ylabel('Diversity')
		# plt.ylabel(key)
		# plt.legend([(1, 'non-recurring'), (2, 'recurring')])
		for pc in parts['bodies']:
			pc.set_facecolor(color)
			pc.set_edgecolor('black')
			pc.set_alpha(1)
		plt.tick_params(
	    axis='x',          # changes apply to the x-axis
	    which='both',      # both major and minor ticks are affected
	    bottom='off',      # ticks along the bottom edge are off
	    top='off',         # ticks along the top edge are off
	    labelbottom='off') # labels along the bottom edge are off
		plt.show()

def _get_stats(master_dir_name, subdir_names):
	stats = {'num_unique_clones': [], 'diversities': [], 'num_bps': []}
	
	for subdir_name in subdir_names:
		full_dir = master_dir_name + subdir_name
		C = np.genfromtxt(full_dir + 'C.tsv', dtype = float)
		U = np.genfromtxt(full_dir + 'U.tsv', dtype = float)
		W = np.genfromtxt(full_dir + 'W.tsv', dtype = int)
		_, l = W.shape
		T = _get_T(open(full_dir + 'T.dot'))

		stats['num_unique_clones'].append(_get_num_unique_clones(T))
		stats['diversities'].append(1.0 - np.sum(U ** 2))
		stats['num_bps'].append(l)

	for key in stats.keys():
		stats[key] = np.array(stats[key])

	return stats

def _get_num_unique_clones(root):
	num_unique_clones = 0
	stack = [root]
	while stack:
		node = stack.pop()
		stack += node.children
		if node.data != '0/0':
			num_unique_clones += 1
	return num_unique_clones

def _get_T(file):
	txt = file.read()
	lines = txt.split('\n')
	lines.pop(0)
	lines.pop(len(lines)-1)
	lines = [ line.replace('\t', '') for line in lines ]
	lines = [ line.replace(' ', '') for line in lines ]
	edges = {} # key is tuple (parent, child). val is edge label
	for line in lines:
		if '->' in line:
			i, j = line.split('->')
			if '[' in j:
				j, edge_label = j.split('[')
				edge_label = '[' + edge_label
				edge_label = edge_label.split('\"')[1]
				edges[(i, j)] = edge_label
	
	nodes = {}
	for (i, j), label in edges.iteritems():
		if i not in nodes.keys():
			nodes[i] = TreeNode(i)
		if j not in nodes.keys():
			nodes[j] = TreeNode(j, data = label)
		nodes[i].children.append(nodes[j])
		nodes[j].parent = nodes[i]

	# find root
	root = nodes[nodes.keys()[0]]
	while root.parent is not None:
		root = root.parent
	return root

def printnow(s, newline = True):
	s = str(s)
	if newline:
		s += '\n'
	sys.stdout.write(s)
	sys.stdout.flush()

#
#   TEMPORARY (REMOVE LATER!!!)
#

def printfull(arr):
	np.set_printoptions(threshold = np.nan)
	print arr


# # # # # # # # # # # # # # # # # # # # # # # # # #
#   C U S T O M   T R E E   N O D E   C L A S S   #
# # # # # # # # # # # # # # # # # # # # # # # # # #

class TreeNode:

	def __init__(self, label, data = ''):
		self.label = label
		self.parent = None
		self.children = []
		self.data = data

	def __str__(self):
		return self.as_str()

	def as_str(self, level = 1):
		indent = '\t'.join([ '' for i in xrange(0, level) ])
		s = indent + 'name:\t' + self.label + '\n'
		s += indent + 'edge:\t' + self.data + '\n'
		s += indent + 'children:\n'
		for child in self.children:
			s += child.as_str(level + 1) + '\n'
		return s


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'template.py', description = "purpose")
	parser.add_argument('-r', '--recur_directory', required = True, type = lambda x: fm.valid_master_dir_with_files_and_ext(parser, x, FNAMES, ''), help = 'directory containing subdirectories each for single patient with recurrence. each subdirectory must contain the files: ' + ', '.join(FNAMES))
	parser.add_argument('-n', '--non_recur_directory', required = True, type = lambda x: fm.valid_master_dir_with_files_and_ext(parser, x, FNAMES, ''), help = 'directory containing subdirectories each for single patient with no recurrence. each subdirectory must contain the files: ' + ', '.join(FNAMES))
	return vars(parser.parse_args(argv))


# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
	main(sys.argv[1:])
