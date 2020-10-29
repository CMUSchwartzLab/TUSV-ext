
#     file: gene_prof.py
#   author: Jingyi Wang
#  created: 10/05/2017
#  modified: 10/09/2017
#  purpose: Genome Prof class. GenomeProf is a profile for each patient that contains 48 chrom profiles (two for each pair of chromsome)


##########
# Import #
##########
import chrm_prof as chpr
import random
import numpy as np
import sys
import os
import copy
import time

sys.path.insert(0, 'helper/')

import combine_copy_nums as ccn

#############
# Functions #
#############

def get_default_mutCount_dict(chrom_dict):
	mutCount_dict = dict()
	for key in chrom_dict:
		mutCount_dict[key] = 0
	return mutCount_dict


def printnow(s, newline = True):
    s = str(s)
    if newline:
        s += '\n'
    sys.stdout.write(s)
    sys.stdout.flush()
##################
# GeneProf class #
##################

class GeneProf:

	# chrom_dict: dicionary
	#             key: a tuple indicates chrom idx and paternal/maternal. (eg. (3, 0): No.3 chromosome from father)
	#             val: a ChrmProf object
	# constants_dict contains: mut_types (list), exp_mut_size (int), exp_mut_count (int/float), cov (int), read_len (int)
	# mutCount: number of mutations for the sample (int)
	# maxCount: total number of mutations for the sample (int)
	# copy_num_dict: dictionary

	def __init__(self, chrom_dict, constants_dict):
		self.chrom_dict = chrom_dict
		self.constants_dict = constants_dict
		self.get_constants()
		self.mutCount = 0
		self.maxCount = self.get_mut_count()
		self.copy_num_dict = self.get_copy_nums_dict()
		self.sv_dict = self.get_sv_read_nums_dict(self.cov, self.read_len)


	def get_constants(self):
		self.mut_types = self.constants_dict['mut_types']
		self.exp_mut_size = self.constants_dict['exp_mut_size']
		self.exp_mut_count = self.constants_dict['exp_mut_count']
		self.cov = self.constants_dict['cov']
		self.read_len = self.constants_dict['read_len']



    # get total number of mutations for the sample randomly based on mutation count distribution
	def get_mut_count(self):
		maxCount = int(round(np.random.poisson(self.exp_mut_count)))
		return maxCount


	# get mutation type, position, size, etc. randomly
	def random_mutation(self):
		mut_type = random.choice(self.mut_types)
		mut_chr = random.choice(list(self.chrom_dict.keys()))
		mut_size = int(round(np.random.exponential(self.constants_dict['exp_mut_size'])))
		while mut_size <= 0:
			mut_size = int(round(np.random.exponential(self.constants_dict['exp_mut_size'])))

		temp = self.chrom_dict[mut_chr].n - mut_size
		while temp <= 0:
			mut_type = random.choice(self.mut_types)
			mut_chr = random.choice(list(self.chrom_dict.keys()))
			mut_size = int(round(np.random.exponential(self.constants_dict['exp_mut_size'])))
			temp = self.chrom_dict[mut_chr].n - mut_size

		mut_bgnPos = random.randint(0, temp)
		mut_endPos = mut_bgnPos + mut_size - 1

		return mut_type, mut_chr, mut_size, mut_bgnPos, mut_endPos



	def is_legal_mutation(self, geneprof_list, mut_type, mut_chr, mut_size, mut_bgnPos, mut_endPos):
		mut_chr_list = [(mut_chr[0], 0), (mut_chr[0], 1)]
		if mut_bgnPos > mut_endPos:
			return False
		if self.chrom_dict[mut_chr_list[0]]._is_splitable(mut_bgnPos, mut_endPos) == False:
			return False
		if self.chrom_dict[mut_chr_list[1]]._is_splitable(mut_bgnPos, mut_endPos) == False:
			return False
		for geneprof in geneprof_list:
			# print 'check', geneprof, 'in geneprof_list, paternal chrom'
			if geneprof.chrom_dict[mut_chr_list[0]]._is_splitable(mut_bgnPos, mut_endPos) == False:
				return False
			# print 'check', geneprof, 'in geneprof_list, maternal chrom'
			if geneprof.chrom_dict[mut_chr_list[1]]._is_splitable(mut_bgnPos, mut_endPos) == False:
				return False
		return True


	# geneProf_list contains list of geneProfs
	def get_legal_random_mutation(self, geneprof_list):
		mut_type, mut_chr, mut_size, mut_bgnPos, mut_endPos = self.random_mutation()

		while self.is_legal_mutation(geneprof_list, mut_type, mut_chr, mut_size, mut_bgnPos, mut_endPos) == False:

			mut_type, mut_chr, mut_size, mut_bgnPos, mut_endPos = self.random_mutation()

		return mut_type, mut_chr, mut_size, mut_bgnPos, mut_endPos


    # make a single mutation
	def mutate(self, geneprof_list):

		mut_type, mut_chr, mut_size, mut_bgnPos, mut_endPos = self.get_legal_random_mutation(geneprof_list)
		# print 'mut_type:', mut_type, 'mut_chr:', mut_chr, 'mut_size:', mut_size, 'mut_bgnPos:', mut_bgnPos, 'mut_endPos:', mut_endPos

		if mut_type == 'amp':
			self.chrom_dict[mut_chr].amp(mut_bgnPos, mut_endPos)

		elif mut_type == 'rem':
			self.chrom_dict[mut_chr].rem(mut_bgnPos, mut_endPos)

		else:
			self.chrom_dict[mut_chr].inv(mut_bgnPos, mut_endPos)

		self.copy_num_dict = self.get_copy_nums_dict()
		self.mutCount += 1
		

    # make multiple mutations 
	def multi_mutations(self, geneprof_list):
		while self.mutCount < self.maxCount:
			self.mutate(geneprof_list)


    # copy_num_dict: dictionary
    #                key: chromosome index (str)
    #                val: tuple of begin_idx_list, end_idx_list, copy_number_list (bgns, ends, cps)
	def get_copy_nums_dict(self):
		result = dict()
		usages = [1,1]
		for (idx, pm) in self.chrom_dict:
			if idx not in result:
				(bgns_p, ends_p, cps_p) = self.chrom_dict[(idx, 0)].get_copy_nums()
				(bgns_m, ends_m, cps_m) = self.chrom_dict[(idx, 1)].get_copy_nums()
				triplets = [[bgns_p, ends_p, cps_p], [bgns_m, ends_m, cps_m]]
				[res_bgns, res_ends, res_cps] = ccn.combine_copy_nums(triplets, usages)
				result[idx] = (res_bgns, res_ends, res_cps)
		
		return result


	# print information of the sample
	def print_info(self):
		l = sorted(self.chrom_dict.keys())
		for (idx,pm) in l:
			print '(', idx, ',', pm, '): ', self.chrom_dict[(idx,pm)].chrm

		print 'copy_num_dict:', self.copy_num_dict
		print 'mutCount:', self.mutCount

		print 'sv_reads_dict:'
		sv_dict = self.get_sv_read_nums_dict(self.cov, self.read_len)
		for i in sorted(sv_dict.keys()):
			print 'chromosome index:', i
			for (pos, isLeft) in sorted(sv_dict[i], key=lambda tup: tup[0]):
				print "(", pos, ",", isLeft, "):", sv_dict[i][(pos, isLeft)]



	def print_chrm_seq(self):
		l = sorted(self.chrom_dict.keys())
		for (idx,pm) in l:
			print '(', idx, ',', pm, '): ', self.chrom_dict[(idx,pm)].chrm


	# input: cov (int), read_len (int)
	# output: sv_dict
	#         key: chromsome index (str)
	#         val: dictionary
	#              key: bp tuple (pos(int), isLeft (bool))
	#              val: dictionary
	#                   key: "mate", val: mate_bp tuple (pos(int), isLeft (bool))
	#                   key: "mated_reads", val: number of reads containing both curr and mated bp (int)
	#                   key: "total_reads", val: number of total reads (int)
	#                   key: "copy_num", val: copy number of current bp (int)
	def get_sv_read_nums_dict(self, cov, read_len):
		result = dict()
		for (idx, pm) in self.chrom_dict:
			if idx not in result:
				temp = dict()
				# paternal chrom
				sv_dict_p = self.chrom_dict[(idx, 0)].get_sv_read_nums(cov, read_len)

				# maternal chrom
				sv_dict_m = self.chrom_dict[(idx, 1)].get_sv_read_nums(cov, read_len)

				# combine paternal and maternal chrom sv dict
				repeated = set()
				for (pos, isLeft) in sv_dict_p:
					if (pos, isLeft) not in sv_dict_m:
						temp[(pos, isLeft)] = sv_dict_p[(pos, isLeft)]
					else:
						repeated.add((pos, isLeft))
						if sv_dict_p[(pos, isLeft)]["mate"] == sv_dict_m[(pos, isLeft)]["mate"]:
							temp[(pos, isLeft)]["mate"] = sv_dict_p[(pos, isLeft)]["mate"]
							temp[(pos, isLeft)]["mated_reads"] = sv_dict_p[(pos, isLeft)]["mated_reads"] + sv_dict_m[(pos, isLeft)]["mated_reads"]
							temp[(pos, isLeft)]["total_reads"] = sv_dict_p[(pos, isLeft)]["total_reads"] + sv_dict_m[(pos, isLeft)]["total_reads"]
							temp[(pos, isLeft)]["copy_num"] = sv_dict_p[(pos, isLeft)]["copy_num"] * 0.5 + sv_dict_m[(pos, isLeft)]["copy_num"] * 0.5
						else:
							print "\n"
							print 'Different mate bp in pair of chromosomes!!!'
							print (pos, isLeft), sv_dict_p[(pos, isLeft)]["mate"], sv_dict_m[(pos, isLeft)]["mate"]
							print "\n"
				for (pos, isLeft) in sv_dict_m:
					if (pos, isLeft) not in repeated:
						temp[(pos, isLeft)] = sv_dict_m[(pos, isLeft)]
				result[idx] = temp
		return result


	def deepcopy(self):
		# deep copy self.chrom_dict
		chrom_dict_new = dict()
		for (idx,pm) in list(self.chrom_dict.keys()):
			chrom_dict_new[(idx,pm)] = self.chrom_dict[(idx,pm)].deepcopy()

		constants_dict_new = copy.deepcopy(self.constants_dict)

		gp = GeneProf(chrom_dict_new, constants_dict_new)

		# copy self.mutCount
		mutCount_new = self.mutCount
		gp.mutCount = mutCount_new

		# copy self.maxCount 
		maxCount_new = self.maxCount
		gp.maxCount = maxCount_new

		# deep copy self.copy_num_dict
		copy_num_dict_new = copy.deepcopy(self.copy_num_dict)
		gp.copy_num_dict = copy_num_dict_new

		# deep copy self.sv_dict
		sv_dict_new = copy.deepcopy(self.sv_dict)
		gp.sv_dict = sv_dict_new

		return gp

