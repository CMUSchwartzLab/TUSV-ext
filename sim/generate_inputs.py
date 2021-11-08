# author: Xuecong Fu
# create the input of CNT-MD format

import pickle
import argparse
import glob
import vcf

def write_cntmd_input(cntmd_dict, output_file):
	num_chroms = len(cntmd_dict[1].keys())
	num_samples = len(cntmd_dict.keys())
	num_segments = []
	for chr, segments in cntmd_dict[1].items():
		num_segments.append(len(segments))
	with open(output_file, 'w') as f:
		print('#PARAMS', file=f)
		print(str(num_chroms) + ' #number of chromosomes', file=f)
		print(str(num_samples) + ' #number of samples', file=f)
		print(*num_segments, sep=' ', file=f, end=' #number of segments for each chromosome\n')
		print('#SAMPLES', file=f)
		for sample, chr_dict in cntmd_dict.items():
			print(str(sample) + ' : ', end='', file=f)
			chr_count = 0
			for chr, segments in chr_dict.items():
				print(*segments, sep=' ', file=f, end='')
				chr_count += 1
				if chr_count < num_chroms:
					print(' | ', end='',file=f)
				else:
					print('', file=f)

def vcf2cntmd_dict(vcfs_dir):
	cntmd_dict = {}
	vcf_files = glob.glob(vcfs_dir + '/*.vcf')
	print(vcf_files)
	for vcf_file in vcf_files:
		sample_id = int(vcf_file.split('.')[0][-1])
		if sample_id not in cntmd_dict.keys():
			cntmd_dict[sample_id] = {}
		vcf_reader = vcf.Reader(open(vcf_file, 'r'))
		while True:
			try:
				record = next(vcf_reader)
				if record.ID[:3] != 'cnv':
					break
				if record.CHROM not in cntmd_dict[sample_id].keys():
					cntmd_dict[sample_id][record.CHROM] = []
				cntmd_dict[sample_id][record.CHROM].append(record.genotype('TUMOR')['CN'][0] + record.genotype('TUMOR')['CN'][1])
			except StopIteration:
				break
	return cntmd_dict

def main1():
	parser = argparse.ArgumentParser()
	parser.add_argument('--cntmd_dir', help='Directory of cntmd input')
	args = parser.parse_args()
	with open(args.cntmd_dir +'/cntmd_dict.pickle', 'rb') as f:
		cntmd_dict = pickle.load(f)
	write_cntmd_input(cntmd_dict, args.cntmd_dir + '/input.samples')

def main2():
	parser = argparse.ArgumentParser()
	parser.add_argument('--vcfs_dir', help='Directory of vcf samples')
	parser.add_argument('--cntmd_dir', help='Directory of cntmd input')
	args = parser.parse_args()
	print(args)
	cntmd_dict = vcf2cntmd_dict(args.vcfs_dir)
	print(cntmd_dict)
	write_cntmd_input(cntmd_dict, args.cntmd_dir + '/input.samples')

if __name__ == '__main__':
	main2()