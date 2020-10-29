#     file: vcf_help.py
#   author: Jesse Eaton
#  created: 12/14/2017
# modified: 12/19/2017
#  purpose: Helper module for writing .vcf files


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import datetime # for changing the created data in .vcf header
import vcf      # install with "pip install pyvcf"


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

BULK_TUMOR_PREFIX = 'BULK_TUMOR' # name of bulk tumor column in output vcf record
# CELL_TYPE_PREFIX = 'CT'

# # # # # # # # # # # # # # # # # # # #
#   V C F   W R I T E R   C L A S S   #
# # # # # # # # # # # # # # # # # # # #

class Writer:

	# input: num_samples (int) number of bulk samples
	#        num_clones (int) number of cell types (number of clones)
	#        metadata_fname (str) file name containing metadata for .vcf output
	def __init__(self, num_samples, num_clones, metadata_fname):
		self.bps = {} # key is id. val is BP
		self.cvs = {} # key is id. val is CV
		self.snames = [ BULK_TUMOR_PREFIX + str(i+1) for i in xrange(0, num_samples) ] # name of each bulk tumor sample column
		self.cnames = [ str(i) for i in xrange(0, num_clones) ] # name of each cell type sample column
		self.metadata_fname = metadata_fname

	#  input: chrm (str) chromosome
	#         pos (int) position on chromosome
	#         ext_left (bool) True if break end extends to the left on original chromosome
	#         rec_id (str) unique identifier for this breakpoint record
	#         mate_id (str) unique identifier of the record for the mate of this breakpoint
	#         mixfs (list of float) mixed copy number of breakpoint in each bulk sample
	#         cps (list of int) [self.num_clones] copy number of breakpoint for each clone
	def add_bp(self, chrm, pos, ext_left, rec_id, mate_id, mixfs, cps):
		self.bps[rec_id] = BP(chrm, pos, ext_left, rec_id, mate_id, mixfs, cps)
		if mate_id in self.bps:
			_set_mates(self.bps[rec_id], self.bps[mate_id])

	#  input: chrm (str) chromosome
	#         bgn (int) beginning position of segment
	#         end (int) ending position of segment
	#         rec_id (str) unique identifier for this copy number segment
	#         mixfs (list of float) mixed copy number of breakpoint in each of the bulk samples
	#         cps (list of int) [self.num_clones] copy number of segment for each clone
	def add_cv(self, chrm, bgn, end, rec_id, mixfs, cps):
		self.cvs[rec_id] = CV(chrm, bgn, end, rec_id, mixfs, cps)

	# input: file (file) file to write all records to
	def write(self, file):
		metadata_file = open(self.metadata_fname, 'r')
		reader = vcf.Reader(metadata_file, 'r')
		reader.samples = self.snames + self.cnames
		reader.metadata['filedate'][0] = datetime.datetime.now().date().strftime('%Y%m%d') # set date to current date
		writer = vcf.Writer(file, reader)

		muts = _merge_dicts(self.bps, self.cvs).values()
		muts.sort(key = lambda mut: mut.pos, reverse = False)
		muts.sort(key = lambda mut: mut.chrm, reverse = False)
		for mut in muts:
			writer.write_record(mut.as_rec(self))
		metadata_file.close()


# # # # # # # # # # # #
#   B P   C L A S S   #
# # # # # # # # # # # #

class BP:

	#  input: chrm (str) chromosome
	#         pos (int) position on chromosome
	#         ext_left (bool) True if break end extends to the left on original chromosome
	#         rec_id (str) unique identifier for this breakpoint record
	#         mate_id (str) unique identifier of the record for the mate of this breakpoint
	#         mixfs (list of float) mixed copy number of breakpoint in each of the bulk samples
	#         cps (list of int) [num_clones] copy number of breakpoint for each clone
	def __init__(self, chrm, pos, ext_left, rec_id, mate_id, mixfs, cps):
		self.chrm = chrm
		self.pos = pos
		self.ext_left = ext_left
		self.rec_id = rec_id
		self.mixfs = mixfs
		self.cps = cps
		self.mate_id = mate_id
		self.mate = None

	#  input: w (Writer)
	# output: rec (vcf._Record)
	def as_rec(self, w):
		mate = self.mate
		ref, qual, filt, fmt = '.', None, [], 'GT:CNADJ'
		alts = [vcf.parser._Breakend(mate.chrm, mate.pos, mate.ext_left, remoteOrientation = self.ext_left, withinMainAssembly = True, connectingSequence = '')]
		info = { 'SVTYPE': 'BND', 'MATEID': mate.rec_id }

		calls = _get_calls(self.mixfs, self.cps, w.snames, w.cnames)
		return vcf.model._Record(self.chrm, self.pos, self.rec_id, ref, alts, qual, filt, info, fmt, w.snames + w.cnames, calls)


# # # # # # # # # # # #
#   C V   C L A S S   #
# # # # # # # # # # # #

class CV:

	#  input: chrm (str) chromosome
	#         bgn (int) beginning position of segment
	#         end (int) ending position of segment
	#         rec_id (str) unique identifier for this copy number segment
	#         mixfs (list of float) mixed copy number of breakpoint in each of the bulk samples
	#         cps (list of int) [self.num_clones] copy number of segment for each clone
	def __init__(self, chrm, bgn, end, rec_id, mixfs, cps):
		self.chrm = chrm
		self.bgn = bgn
		self.end = end
		self.rec_id = rec_id
		self.mixfs = mixfs
		self.cps = cps
		self.pos = bgn # for integrated sorting with SVs

	#  input: w (Writer)
	# output: rec (vcf._Record)
	def as_rec(self, w):
		ref, qual, filt, fmt = '.', None, [], 'GT:CN'
		alts = [vcf.model._SV('CNV')]
		info = { 'IMPRECISE': True, 'END': self.end }
		calls = _get_calls(self.mixfs, self.cps, w.snames, w.cnames)
		return vcf.model._Record(self.chrm, self.pos, self.rec_id, ref, alts, qual, filt, info, fmt, w.snames + w.cnames, calls)


# # # # # # # # # # # # # # # # # #
#   H E L P E R   C L A S S E S   #
# # # # # # # # # # # # # # # # # #

class BpCallData:
	def __init__(self, cnadj, gt = '1|0'):
		self.GT = gt
		self.CNADJ = str(cnadj)
		self.__getitem__ = self

	def __call__(self, var):
		return [self.CNADJ]


# # # # # # # # # # # # # # # # # # # # #
#   P R I V A T E   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # #

def _set_mates(bp1, bp2):
	bp1.mate = bp2
	bp2.mate = bp1

def _merge_dicts(d1, d2):
	d3 = {}
	for k, v in d1.iteritems():
		d3[k] = v
	for k, v in d2.iteritems():
		d3[k] = v
	return d3

#  input: mixfs (list of float) [num_bulk_samples] mixed copy number of breakpoint in each of the bulk samples
#         cps (list of int) [num_clones] copy number of segment for each clone
#         snames (list of str) [num_bulk_samples] name of each bulk tumor sample column
#         cnames (list of str) [num_clones] name of each cell type sample column
def _get_calls(mixfs, cps, snames, cnames):
	calls = []
	for i, f in enumerate(mixfs):
		calls.append(vcf.model._Call(i, snames[i], BpCallData(cnadj = f)))
	for i, cp in enumerate(cps):
		calls.append(vcf.model._Call(i, cnames[i], BpCallData(cnadj = int(cp))))
	return calls
