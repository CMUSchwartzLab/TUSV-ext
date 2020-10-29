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

from os.path import isfile, join


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

# creates file if file does not already exist
def touch(fname, times = None):
	with open(fname, 'a'):
		os.utime(fname, times)

def is_valid_file(parser, arg):
	if not os.path.exists(arg):
		parser.error('The file \"' + str(arg) + '\" could not be found.')
	else:
		return str(arg)

def append_to_file(fname, msg):
	touch(fname)
	with open(fname, "a") as f:
		f.write(msg)

def mkdir(dname):
	if not os.path.exists(dname):
		os.makedirs(dname)

# returns string as directory. adds error to parser if no valid directory
def valid_dir(parser, arg):
	if not os.path.exists(arg):
		parser.error('The directory \"' + str(arg) + '\" could not be found.')
	return _directorize(arg)

#  input: parser (argparser.parser)
#         arg (str) full path of directory
#         ext (str) extension (ex: .vcf). directory must have at least one of these files
# output: arg (str) full path of directory with '/' as needed
def valid_dir_ext(parser, arg, ext):
	if not os.path.exists(arg):
		parser.error('The directory \"' + str(arg) + '\" could not be found.')
	arg = _directorize(arg)
	if len(_fnames_with_extension(arg, ext)) == 0:
		parser.error('The directory \"' + str(arg) + '\" contained no ' + str(ext) + ' files.')
	return arg

#  input: parser (argparser.parser)
#         arg (str) full path of directory
#         fnames (list of str) names of files in directory. directory must contain all these files
# output: arg (str) full path of directory with '/' as needed
def valid_dir_with_files(parser, arg, fnames):
	if not os.path.exists(arg):
		parser.error('The directory \"' + str(arg) + '\" could not be found.')
	arg = _directorize(arg)
	fnames_not_found = _fnames_not_found(arg, fnames)
	if len(fnames_not_found) > 0:
		parser.error('The directory \"' + str(arg) + '\" did not contain ' + ', '.join(fnames_not_found) + ' file(s).')
	return arg

#  input: parser (argparser.parser)
#         arg (str) full path of directory contining subdirectories
#         fnames (list of str) names of files in subdirectories. subdirectories must contain all these files
#         ext (str) extension (ex: .vcf). subdirectories must have at least one of these files
# output: arg (str) full path of directory with '/' as needed
def valid_master_dir_with_files_and_ext(parser, arg, fnames, ext):
	if not os.path.exists(arg):
		parser.error('The directory \"' + str(arg) + '\" could not be found.')
	arg = _directorize(arg)
	for subdir, dirs, _ in _walklevel(arg):
		if not dirs:
			parser.error('The directory \"' + str(subdir) + '\" has no subdirectories.')
		for d in dirs:
			valid_dir_with_files(parser, arg + d, fnames)
			valid_dir_ext(parser, arg + d, ext)
	return arg


#
#   File creation functions
#

# copies file structure (only one level deep) from in_dir to out_dir
def cp_file_structure_to_out_dir(in_dir, out_dir):
	for _, subdirs, _ in _walklevel(in_dir):
		ot_subdirs = [ _directorize(out_dir + d) for d in subdirs ]
		for ot_subdir in ot_subdirs:
			if not os.path.exists(ot_subdir):
				os.makedirs(ot_subdir)

def get_subdir_names(d):
	for _, subdirs, _ in _walklevel(d):
		return [ _directorize(s) for s in subdirs ]

def get_fnames_in_dir(d):
	return [f for f in os.listdir(d) if isfile(join(d, f))]

def get_fnames_in_dir_with_ext(d, ext):
	fnames = get_fnames_in_dir(d)
	out = []
	for fname in fnames:
		if fname.endswith(ext):
			out.append(fname)
	return out

#
#   non-file maninging input functions
#

def valid_int_in_range(parser, arg, lo, hi):
	try:
		arg = int(arg)
	except:
		parser.error(str(arg) + ' must be an integer')
	if arg < lo or arg > hi:
		parser.error(str(arg) + ' must be an integer between ' + str(lo) + ' and ' + str(hi))
	return arg

def valid_float_above(parser, arg, lo):
	try:
		arg = float(arg)
	except:
		parser.error(str(arg) + ' must be a float')
	if arg < lo:
		parser.error(str(arg) + ' must be a float above ' + str(lo))
	return arg


# # # # # # # # # # # # # # # # # # # #
#   H E L P E R   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # #

# add "/" to end of directory name if necessary
def _directorize(dir_name):
	if dir_name.endswith('/'):
		return dir_name
	return dir_name + '/'

# os.walk but only goes recursively down level levels
def _walklevel(some_dir, level = 0):
	some_dir = some_dir.rstrip(os.path.sep)
	assert os.path.isdir(some_dir)
	num_sep = some_dir.count(os.path.sep)
	for root, dirs, files in os.walk(some_dir):
		yield root, dirs, files
		num_sep_this = root.count(os.path.sep)
		if num_sep + level <= num_sep_this:
			del dirs[:]

# returns all files in the directory with extension ext (ex. ext = '.vcf')
def _fnames_with_extension(directory, ext):
	files = []
	for file in os.listdir(directory):
		if file.endswith(ext):
			files.append(file)
	return files

def _fnames_not_found(directory, fnames):
	fnames_found = [ fname for fname in os.listdir(directory) ]
	return list(set(fnames) - set(fnames_found))
