#     file: multi_tusv.py
#   author: Jesse Eaton
#  created: 12/3/2017
# modified: 12/3/2017
#  purpose: Runs unmixing mixed copy numbers for breakpoints and segments and infers phylogeny
#             with phylogenetic constraints across multiple patients


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments

# custom modules
import tusv
sys.path.insert(0, 'help/')
import file_manager as fm
import printer as pt


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)

	in_dir, out_dir = args['input_directory'], args['output_directory']
	fm.cp_file_structure_to_out_dir(in_dir, out_dir)
	for subdir_name in sorted(fm.get_subdir_names(in_dir)):
		sub_in_dir = in_dir + subdir_name
		sub_out_dir = out_dir + subdir_name
		pt.printnow(''.join([ '\n' for _ in xrange(0, 10) ]))
		pt.printnow(' '.join([ '=' for _ in xrange(0, 30) ]))
		msg = 'Running: ' + subdir_name if not os.listdir(sub_out_dir) else 'ALREADY RAN: ' + subdir_name
		pt.printnow('\t' + msg)
		pt.printnow(' '.join([ '=' for _ in xrange(0, 30) ]))

		if not os.listdir(sub_out_dir): # directory is empty
			tusv.unmix(sub_in_dir, sub_out_dir, args['num_leaves'], args['c_max'], args['lambda1'], args['lambda2'], args['restart_iters'], args['cord_desc_iters'], args['processors'], args['time_limit'], args['metadata_file'], args['num_subsamples'], args['overide_lambdas'])

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'multi_tusv.py', description = "runs tusv.py multiple times. Unmixes mixed copy numbers for breakpoints and segments and infers phylogenies with phylogenetic constraints across multiple patients")
	parser.add_argument('-i', '--input_directory', required = True, type = lambda x: fm.valid_master_dir_with_files_and_ext(parser, x, [], '.vcf'), help = 'directory containing multiple subdirectories each containing a .vcf for each sample from a single patient')
	parser.add_argument('-o', '--output_directory', required = True, type = lambda x: fm.valid_dir(parser, x), help = 'empty directory where multiple subdirectories will be created. each will then contain output U.tsv, C.tsv, and T.dot files')
	tusv.set_non_dir_args(parser)
	return vars(parser.parse_args(argv))

# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
	main(sys.argv[1:])
