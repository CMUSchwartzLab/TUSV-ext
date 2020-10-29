
# INPUT: 2017_09_18_metadata.vcf 
# output: three samples(.vcf) for patient 1 in dummy folder
# Author: Jingyi Wang
# Date: 2017_09_18

##################
##### IMPORT #####
##################

import sys      
import os       
import argparse 
import vcf      
import operator
import datetime
import shutil

#####################
##### FUNCTIONS #####
#####################

def main(argv):
	args = get_args(argv)
	metaFile = args['inputFile']
	vcf_reader = vcf.Reader(open(metaFile, 'r'))
	# print dir(vcf_reader.metadata)
	# print vcf_reader.metadata
	vcf_reader.metadata['filedate'][0] = datetime.datetime.now().date().strftime('%Y%m%d') # set date to current date

	# create patient 1 folder
	script_dir = os.path.dirname(os.path.abspath(__file__))
	outputFolder = script_dir + '/dummy' + '/patient1'
	if os.path.exists(outputFolder):
		shutil.rmtree(outputFolder)
	os.makedirs(outputFolder)
	
	# hard code three samples for patient 1
	# sample 1
	file1 = outputFolder + '/sample1.vcf'
	vcf_writer1 = vcf.Writer(open(file1, 'w'), vcf_reader)
	vcf_writer1.write_record(generate_sv('1', 20000001, 'sv0', '1', 50000000, True, False, str(), True, 'sv1', '1|0', '1', 38, 140))
	vcf_writer1.write_record(generate_sv('1', 50000000, 'sv1', '1', 20000001, True, True, str(), True, 'sv0', '1|0', '1', 38, 130))
	vcf_writer1.write_record(generate_cnv('1', 1, 'cnv0', 'CNV', 20000000, '1|1', [1, 1]))
	vcf_writer1.write_record(generate_cnv('1', 20000001, 'cnv1', 'CNV', 50000000, '1|1', [2, 1]))
	vcf_writer1.write_record(generate_cnv('1', 50000001, 'cnv2', 'CNV', 100000000, '1|1', [1, 1]))
	# vcf_writer1.write_record(generate_cnv('2', 10000001, 'cnv1', 'CNV', 20000000, '1|1', [3, 1]))
	# vcf_writer1.write_record(generate_cnv('3', 1000, 'cnv2', 'CNV', 1500, '1|1', [3, 0]))

	# sample 2
	file2 = outputFolder + '/sample2.vcf'
	vcf_writer2 = vcf.Writer(open(file2, 'w'), vcf_reader)
	vcf_writer2.write_record(generate_sv('1', 20000000, 'sv0', '1', 50000000, True, False, str(), True, 'sv2', '1|0', '1', 13, 97))
	vcf_writer2.write_record(generate_sv('1', 50000000, 'sv2', '1', 20000000, True, False, str(), True, 'sv0', '1|0', '1', 13, 113))
	vcf_writer2.write_record(generate_sv('1', 20000001, 'sv1', '1', 50000001, True, True, str(), True, 'sv3', '1|0', '1', 17, 98))
	vcf_writer2.write_record(generate_sv('1', 50000001, 'sv3', '1', 20000001, True, True, str(), True, 'sv1', '1|0', '1', 17, 104))
	vcf_writer2.write_record(generate_cnv('1', 1, 'cnv0', 'CNV', 20000000, '1|1', [1, 1]))
	vcf_writer2.write_record(generate_cnv('1', 20000001, 'cnv1', 'CNV', 40000000, '1|1', [3, 1]))
	vcf_writer2.write_record(generate_cnv('1', 40000001, 'cnv2', 'CNV', 50000000, '1|1', [1, 1]))
	vcf_writer2.write_record(generate_cnv('1', 50000001, 'cnv3', 'CNV', 100000000, '1|1', [1, 1]))

	# sample 3
	file3 = outputFolder + '/sample3.vcf'
	vcf_writer3 = vcf.Writer(open(file3, 'w'), vcf_reader)
	vcf_writer3.write_record(generate_sv('1', 20000001, 'sv0', '1', 50000000, True, False, str(), True, 'sv1', '1|0', '1', 21, 100))
	vcf_writer3.write_record(generate_sv('1', 50000000, 'sv1', '1', 20000001, True, True, str(), True, 'sv0', '1|0', '1', 21, 101))
	vcf_writer3.write_record(generate_sv('1', 60000000, 'sv2', '1', 80000001, True, True, str(), True, 'sv3', '1|0', '1', 37, 167))
	vcf_writer3.write_record(generate_sv('1', 80000001, 'sv3', '1', 60000000, True, False, str(), True, 'sv2', '1|0', '1', 37, 173))
	vcf_writer3.write_record(generate_cnv('1', 1, 'cnv0', 'CNV', 20000000, '1|1', [1, 1]))
	vcf_writer3.write_record(generate_cnv('1', 20000001, 'cnv1', 'CNV', 50000000, '1|1', [2, 1]))
	vcf_writer3.write_record(generate_cnv('1', 50000001, 'cnv2', 'CNV', 60000000, '1|1', [1, 1]))
	vcf_writer3.write_record(generate_cnv('1', 60000001, 'cnv3', 'CNV', 80000000, '1|1', [1, 0]))
	vcf_writer3.write_record(generate_cnv('1', 80000001, 'cnv4', 'CNV', 100000000, '1|1', [1, 1]))


# type(chrom): str
# type(pos): int
# type(rec_id): str
# ref = str()
# type(alts): list
#              type(alts[0]) = class 'vcf.model._Breakend'
#              dir(alts[0]): [..., 'chr', 'connectingSequence', 'orientation', 'pos', 'remoteOrientation', 'type', 'withinMainAssembly']
#                           eg. alts[0] = ]1:149965077]
#                               'chr' = str(1), 'connectingSequence' = str(), 'orientation' = True, 'pos' = int(149965077)
#                               'remoteOrientation' = False, 'type' = 'BND', 'withinMainAssembly' = True
# qual = None
# filter = list()
# type(info): dict, info['SVTYPE'] = 'BND', info['MATEID'] = mate_sv_rec_id (str)
# fmt = 'GT:CNADJ'
# sample = ['TUMOR', 'NORMAL']
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
	calls = [vcf.model._Call(0, 'TUMOR', svCallData(gt,cnadj,bdp,dp)), vcf.model._Call(1, 'NORMAL', svCallData('0|0','0', 0, dp))]
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


###################
##### CLASSES #####
###################

class svCallData:
	def __init__(self, gt = '0|0', cnadj = '0', bdp = '0', dp = '100'):
		self.GT = gt
		self.CNADJ = cnadj
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
##### COMMAND LINE ARGUMENT FUNCTION #####
###########################################

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'genetate_dummy.py', description = "generate one dummy patient with three samples")
	parser.add_argument('-f', '--input_file', type = str, dest = "inputFile", required = True)
	# parser.add_argument('-o', '--output_folder', type = str, dest = "outputFolder", required = True)
	return vars(parser.parse_args(argv))


##############################
##### CALL MAIN FUNCTION #####
##############################

if __name__ == "__main__":
	main(sys.argv[1:])
