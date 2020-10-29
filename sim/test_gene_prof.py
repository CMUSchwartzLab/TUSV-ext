
#     file: test_gene_prof.py
#   author: Jingyi Wang
#  created: 10/09/2017
#  modified: 10/12/2017
#  purpose: test gene_prof.py 


##########
# Import #
##########

import gene_prof as gnpr
import chrm_prof as chpr


#############
# Test case 1
#############

chrom_dict = dict()
chrom_dict[('1', 0)] = chpr.ChrmProf("AAABBBCCCDDDEEEFFF")
chrom_dict[('1', 1)] = chpr.ChrmProf("AAABBBCCCDDDEEEFFF")

chrom_dict[('2', 0)] = chpr.ChrmProf("RRRRRSSSSSTTTTTUUUUUVVVVV")
chrom_dict[('2', 1)] = chpr.ChrmProf("RRRRRSSSSSTTTTTUUUUUVVVVV")

chrom_dict[('3', 0)] = chpr.ChrmProf("XXXXYYYYZZZZ")
chrom_dict[('3', 1)] = chpr.ChrmProf("XXXXYYYYZZZZ")

constants_dict = dict()
constants_dict['mut_types'] = ['amp', 'rem', 'inv']
constants_dict['mut_size_mean'] = 5
constants_dict['mut_size_var'] = 2
constants_dict['mut_count_mean'] = 3
constants_dict['mut_count_var'] = 1
constants_dict['cov'] = 20
constants_dict['read_len'] = 5



gp = gnpr.GeneProf(chrom_dict, constants_dict)
gp.print_info()
print ''

<<<<<<< HEAD
gp.multi_mutations()
=======
gp.muti_mutations()
>>>>>>> c5222f8eb70a1e1bff9819725c4f5f8e766f2e6e
print '\n\n'
print 'gp info:'
gp.print_info()

print '\n\n'
gp_copied = gp.deepcopy()

gp_copied.mutate()
print 'gp copied info:'
gp_copied.print_info()

print '\n\n'
print 'gp info:'
gp.print_info()


