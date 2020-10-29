#     file: experiment.py
#   author: Jesse Eaton
#  created: 10/24/2017
# modified: 10/24/2017
#  purpose: runs tusv.py on mulitple patients and validates the results


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys        # for command line arguments
import os         # for manipulating files and folders
import argparse   # for command line arguments
import subprocess # for calling tusv.py and validate.py
import collections
import numpy as np
import multiprocessing as mp

sys.path.insert(0, 'help/')
import file_manager as fm   # sanitizes file and directory arguments
import printer as pt

sys.path.insert(0, 'validate/')
import validate as vd

import tusv as tusv


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

FNAMES = ['C.tsv', 'U.tsv', 'T.dot'] # these files must be in 'actual' and 'expected' dirs. MUST BE IN THIS ORDER
EXTENSION = '.vcf'
MAX_NUM_LEAVES = 10
MAX_COPY_NUM = 20
MAX_CORD_DESC_ITERS = 1000
MAX_RESTART_ITERS = 1000
NUM_CORES = mp.cpu_count()
STR_DTYPE = '|S30'


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	fm.cp_file_structure_to_out_dir(args['input_directory'], args['output_directory'])
	subdir_names = fm.get_subdir_names(args['input_directory'])
	tusv.write_readme(args['output_directory'], args, os.path.basename(__file__))

	CBs, CSs, Cs, Us, Ts, FUCs, objs, ls, rs = [], [], [], [], [], [], [], [], [] # scores for copy number of breakpoints, segments, usages and phylogeny
	for subdir_name in subdir_names:
		in_dir = args['input_directory'] + subdir_name
		out_dir = args['output_directory'] + subdir_name

		if not os.listdir(out_dir): # empty directory
			pt.printnow('#\n' * 5 + '\nrunning ' + subdir_name + '\n\n' + '#\n' * 5 + '\n')
			tusv.unmix(in_dir, out_dir, args['num_leaves'], args['c_max'], args['lambda1'], args['lambda2'], args['restart_iters'], args['cord_desc_iters'], args['processors'], args['time_limit'], args['metadata_file'], args['num_subsamples'], args['overide_lambdas'])
		else:
			pt.printnow('#\n' * 5 + '\n\nALREADY RAN ' + subdir_name + '\n\n' + '#\n' * 5 + '\n')

		score_Cb, score_Cs, score_C, score_U, dist_T, score_FUC, obj_val, l, r = vd.get_scores(out_dir, in_dir)
		pt.printnow(' Cb: ' + str(score_Cb))
		pt.printnow(' Cs: ' + str(score_Cs))
		pt.printnow('  C: ' + str(score_C))
		pt.printnow('  U: ' + str(score_U))
		pt.printnow('  T: ' + str(dist_T))
		pt.printnow('FUC: ' + str(score_FUC))
		pt.printnow('obj: ' + str(obj_val))

		CBs.append(score_Cb)
		CSs.append(score_Cs)
		Cs.append(score_C)
		Us.append(score_U)
		Ts.append(dist_T)
		FUCs.append(score_FUC)
		objs.append(obj_val)
		ls.append(l)
		rs.append(r)

	report(args['input_directory'], args['output_directory'], CBs, CSs, Cs, Us, Ts, FUCs, objs, ls, rs, subdir_names)

def report(in_dir, out_dir, CBs, CSs, Cs, Us, Ts, FUCs, obj_vals, ls, rs, subdir_names):
	patient_names = [ os.path.basename(os.path.normpath(dname)) for dname in subdir_names ]
	out_dname = out_dir + 'validation_results/'
	fm.mkdir(out_dname)
	
	# put all scores in .tsv 
	fname = out_dname + 'scores.tsv'
	fm.touch(fname)
	X = np.array([CBs, CSs, Cs, Us, Ts, FUCs, obj_vals, ls, rs], dtype = float).T.astype(STR_DTYPE)
	row_header = np.array(['C_bp_score', 'C_sg_score', 'C_tot_score', 'U_score', 'T_score', 'F_score', 'obj_val', 'l', 'r']).astype(STR_DTYPE)
	col_header = np.array([''] + patient_names).astype(STR_DTYPE)
	X = np.column_stack((col_header, np.ma.row_stack((row_header, X))))
	np.savetxt(fname, X, delimiter = '\t', fmt = '%s')


	fname = out_dname + 'report.txt'
	fm.touch(fname)

	orig_stdout = sys.stdout
	f = open(fname, 'w')
	sys.stdout = f

	pt.printnow(' input directory: ' + in_dir)
	pt.printnow('output directory: ' + out_dir + '\n')

	results = [
		('|C_tru - C_obs| breakpoints only', CBs),
		('|C_tru - C_obs| segments only', CSs),
		('|C_tru - C_obs|', Cs),
		('|U_tru - U_obs|', Us),
		('|T_tru - T_obs| robinson foulds distane', Ts),
		('|F - UC|', FUCs),
		('|F - UC| + L1*tree_cost + L2*', obj_vals)
	]

	for msg, vals in results:
		pt.printnow(msg)
		pt.printnow('  all: ' + ', '.join([ str(val) for val in sorted(vals) ]))
		pt.printnow('  avg: ' + str(np.average(vals)))
		pt.printnow('  std: ' + str(np.std(vals)))
		pt.printnow('')

	sys.stdout = orig_stdout
	f.close()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'experiment.py', description = "runs tusv.py on mulitple patients and validates the results")
	parser.add_argument('-i', '--input_directory', required = True, type = lambda x: fm.valid_master_dir_with_files_and_ext(parser, x, FNAMES, EXTENSION), help = 'directory one or multiple subdirectories. each subdirectory should contain T.dot, C.tsv, and U.tsv files that were generated by sim.py')
	parser.add_argument('-o', '--output_directory', required = True, type = lambda x: fm.valid_dir(parser, x), help = 'directory where results from each experiment will go')
	tusv.set_non_dir_args(parser)
	return vars(parser.parse_args(argv))

def is_valid_file(parser, arg):
	if not os.path.exists(arg):
		parser.error('The file \"' + str(arg) + '\" could not be found.')
	else:
		return open(arg, 'r')


# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
	main(sys.argv[1:])
