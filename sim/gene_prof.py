
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
import operator

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

	def random_mutation_snv(self):
		mut_num = np.random.poisson(self.constants_dict["snv_mut_lambda"])
		mut_chr = random.choice(list(self.chrom_dict.keys()), size=mut_num)
		mut_pos = np.random.randint(0, self.chrom_dict[mut_chr].n, size=mut_num)
		return mut_num, mut_chr, mut_pos

	def is_legal_trans(self, chr1, ins_Pos, chr2, bgn, end):
		if chr1 == chr2 and bgn <= ins_Pos <= end:
			return False
		if self.chrom_dict[(chr2[0], chr2[1])]._is_splitable_one(ins_Pos) is False:
			return False
		return True

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


    # make a single mutation ###xf: lack translocation
	def mutate(self, geneprof_list, snv):

		mut_type, mut_chr, mut_size, mut_bgnPos, mut_endPos = self.get_legal_random_mutation(geneprof_list)
		# print 'mut_type:', mut_type, 'mut_chr:', mut_chr, 'mut_size:', mut_size, 'mut_bgnPos:', mut_bgnPos, 'mut_endPos:', mut_endPos

		if mut_type == 'amp':
			amp_num = np.random.randint(1, 8)
			self.chrom_dict[mut_chr].amp(mut_bgnPos, mut_endPos, amp_num, snv)

		elif mut_type == 'rem':
			self.chrom_dict[mut_chr].rem(mut_bgnPos, mut_endPos, snv)

		elif mut_type == 'inv':
			self.chrom_dict[mut_chr].inv(mut_bgnPos, mut_endPos, snv)

		elif mut_type == 'trans':
			mut_chr2 = random.choice(list(self.chrom_dict.keys()))
			ins_Pos = random.randint(0, self.chrom_dict[mut_chr2].n)
			while not self.is_legal_trans(mut_chr, ins_Pos, mut_chr2, mut_bgnPos, mut_endPos):
				mut_chr2 = random.choice(list(self.chrom_dict.keys()))
				ins_Pos = random.randint(0, self.chrom_dict[mut_chr2].n)
			self.chrom_dict[mut_chr] = self.chrom_dict[mut_chr2].trans(self.chrom_dict[mut_chr], ins_Pos, mut_bgnPos, mut_endPos, snv)

		self.copy_num_dict = self.get_copy_nums_dict()
		self.mutCount += 1

    # make multiple mutations 
	def multi_mutations(self, geneprof_list):
		if self.constants_dict['snv_mut_lambda'] is not None:
			snv=True
		else:
			snv=False
		while self.mutCount < self.maxCount:
			self.mutate(geneprof_list, snv)
		if self.constants_dict['snv_mut_lambda'] is not None:
			mut_num_snv, mut_chr_snv, mut_pos_snv = self.random_mutation_snv()
			for i in range(mut_num_snv):
				self.chrom_dict[mut_chr_snv[i]].point_mutation(mut_pos_snv[i])


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
				[res_bgns, res_ends, res_cps_1, res_cps_2] = ccn.combine_copy_nums_allelic(triplets, usages) ### xf: place where allelic CNVs are combined, modified the function to be allelic specific CNs
				result[idx] = (res_bgns, res_ends, res_cps_1, res_cps_2)
		
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

	### xf: add get snv dict
	def get_snv_dict(self):
		snvs = dict()
		for (idx, pm) in self.chrom_dict:
			snvs = self.chrom_dict[(idx, pm)].get_snvs(snvs)
		return snvs

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
	def get_sv_read_nums_dict(self, cov, read_len):  ### xf: generate sv for both alleles separately and then combine them
		result = dict()
		others = dict()
		for (idx, pm) in self.chrom_dict:
			if pm == 1:
				continue
			if idx not in result:
				temp = dict()
				# paternal chrom
				sv_dict_p, others_p = self.chrom_dict[(idx, 0)].get_sv_read_nums(cov, read_len, idx, 0)
				for tup_chr, items in others_p.items():
					if items == {}:
						continue
					if tup_chr not in others.keys():
						others[tup_chr] = items
					else:
						for (tup_pos, value) in items.items():
							if tup_pos not in others[tup_chr].keys():
								others[tup_chr][tup_pos] = value
							else:
								print(others[tup_chr][tup_pos]["mate"] == others_p[tup_chr][tup_pos]["mate"])
								others[tup_chr][tup_pos]["mated_reads"] += others_p[tup_chr][tup_pos]["mated_reads"]
								others[tup_chr][tup_pos]["total_reads"] += others_p[tup_chr][tup_pos]["total_reads"]
								others[tup_chr][tup_pos]["copy_num"] += others_p[tup_chr][tup_pos]["copy_num"]


				# maternal chrom
				sv_dict_m, others_m = self.chrom_dict[(idx, 1)].get_sv_read_nums(cov, read_len, idx, 1)
				for tup_chr, items in others_m.items():
					if items == {}:
						continue
					if tup_chr not in others.keys():
						others[tup_chr] = items
					else:
						for (tup_pos, value) in items.items():
							if tup_pos not in others[tup_chr].keys():
								others[tup_chr][tup_pos] = value
							else:
								print(others[tup_chr][tup_pos]["mate"] == others_m[tup_chr][tup_pos]["mate"])
								others[tup_chr][tup_pos]["mated_reads"] += others_m[tup_chr][tup_pos]["mated_reads"]
								others[tup_chr][tup_pos]["total_reads"] += others_m[tup_chr][tup_pos]["total_reads"]
								others[tup_chr][tup_pos]["copy_num"] += others_m[tup_chr][tup_pos]["copy_num"]
				# combine paternal and maternal chrom sv dict
				repeated = set()
				for (pos, isLeft, chr_, pm_) in sv_dict_p:
					if (pos, isLeft, chr_, 1) not in sv_dict_m:
						temp[(pos, isLeft, chr_)] = sv_dict_p[(pos, isLeft, chr_, 0)]
						temp[(pos, isLeft, chr_)]["mate"] = (
						sv_dict_p[(pos, isLeft, chr_, 0)]["mate"][0], sv_dict_p[(pos, isLeft, chr_, 0)]["mate"][1],
						sv_dict_p[(pos, isLeft, chr_, 0)]["mate"][2])
					else:
						repeated.add((pos, isLeft, chr_))
						if sv_dict_p[(pos, isLeft, chr_, 0)]["mate"][:3] == sv_dict_m[(pos, isLeft, chr_, 1)]["mate"][:3]:
							temp[(pos, isLeft, chr_)]["mate"] = (sv_dict_p[(pos, isLeft, chr_, 0)]["mate"][0], sv_dict_p[(pos, isLeft, chr_, 0)]["mate"][1],
																 sv_dict_p[(pos, isLeft, chr_, 0)]["mate"][2])
							temp[(pos, isLeft, chr_)]["mated_reads"] = sv_dict_p[(pos, isLeft, chr_, 0)]["mated_reads"] + sv_dict_m[(pos, isLeft, chr_, 1)]["mated_reads"]
							temp[(pos, isLeft, chr_)]["total_reads"] = sv_dict_p[(pos, isLeft, chr_, 0)]["total_reads"] + sv_dict_m[(pos, isLeft, chr_, 1)]["total_reads"]
							temp[(pos, isLeft, chr_)]["copy_num"] = sv_dict_p[(pos, isLeft, chr_, 0)]["copy_num"] + sv_dict_m[(pos, isLeft, chr_, 1)]["copy_num"]
							###xf: change the relative copy number (?) to total copy number
						else:
							print "\n"
							print 'Different mate bp in pair of chromosomes!!!'
							print (pos, isLeft, chr_, pm_), sv_dict_p[(pos, isLeft, chr_, 0)]["mate"], sv_dict_m[(pos, isLeft, chr_, 1)]["mate"]
							print "\n"
				for (pos, isLeft, chr_, pm_) in sv_dict_m:
					if (pos, isLeft, chr_) not in repeated:
						temp[(pos, isLeft, chr_)] = sv_dict_m[(pos, isLeft, chr_, 1)]
						temp[(pos, isLeft, chr_)]["mate"] = (sv_dict_m[(pos, isLeft, chr_, 1)]["mate"][0], sv_dict_m[(pos, isLeft, chr_, 1)]["mate"][1],
															 sv_dict_m[(pos, isLeft, chr_, 1)]["mate"][2])
				result[idx] = temp
		for (chr, pm) in others.keys():
			for (pos, isLeft, chr_, pm_) in others[(chr, pm)]:
				print (pos, isLeft, chr_, pm_)
				assert chr_ == chr
				assert pm_ == pm
				if (pos, isLeft, chr) in result[chr].keys():
					if result[chr][(pos, isLeft, chr)]["mate"] == others[(chr, pm)][(pos, isLeft, chr_, pm_)]["mate"][0:3]:
						result[chr][(pos, isLeft, chr)]["mated_reads"] += others[(chr, pm)][(pos, isLeft, chr_, pm_)]["mated_reads"]
						result[chr][(pos, isLeft, chr)]["total_reads"] += others[(chr, pm)][(pos, isLeft, chr_, pm_)]["total_reads"]
						result[chr][(pos, isLeft, chr)]["copy_num"] += others[(chr, pm)][(pos, isLeft, chr_, pm_)]["copy_num"]
					else:
						print("not consist",result[chr][(pos, isLeft, chr)]["mate"], others[(chr, pm)][(pos, isLeft, chr_, pm_)]["mate"][0:3])
				else:
					result[chr][(pos, isLeft, chr)] = others[(chr, pm)][(pos, isLeft, chr_, pm_)]
					result[chr][(pos, isLeft, chr)]["mate"] = (others[(chr, pm)][(pos, isLeft, chr_, pm_)]["mate"][0],
						others[(chr, pm)][(pos, isLeft, chr_, pm_)]["mate"][1], others[(chr, pm)][(pos, isLeft, chr_, pm_)]["mate"][2])
		return result


	def deepcopy(self):
		# deep copy self.chrom_dict
		chrom_dict_new = dict()
		other_muts_dict = {}
		muts_dict = {}
		for (idx,pm) in list(self.chrom_dict.keys()):
			chrom_dict_new[(idx,pm)], muts, other_muts_dict = self.chrom_dict[(idx,pm)].deepcopy_(other_muts_dict)
			muts_dict[(idx,pm)] = muts
		for (idx,pm) in list(self.chrom_dict.keys()):
			if (idx,pm) in other_muts_dict.keys():
				print("others",other_muts_dict[(idx,pm)])
				muts_dict[(idx,pm)] += other_muts_dict[(idx,pm)]
			muts_sorted = sorted(muts_dict[(idx,pm)], key = lambda x: x.bgn)
			n = len(muts_sorted)
			for i in xrange(0, n):
				if i != 0:
					muts_sorted[i].l = muts_sorted[i-1]
				if i != n-1:
					muts_sorted[i].r = muts_sorted[i+1]
			chrom_dict_new[(idx,pm)].mut = muts_sorted[0]

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

