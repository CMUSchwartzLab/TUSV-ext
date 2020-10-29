
#     file: test_sim.py
#   author: Jingyi Wang
#  created: 10/17/2017
#  modified: 10/12/2017
#  purpose: test sim.py 



##########
# Import #
##########

import gene_prof as gnpr
import chrm_prof as chpr
import sim as sm


####################
# Test random tree #
####################

n = 5
repeat_time = 3

# given n (int), generate repeat_time(int) random trees and their corresponding edge lists
def test_random_tree(n, repeat_time):

	print 'Test random tree:'
	for i in range(repeat_time):
		print 'test', i
		l = sm.random_get_tree(n)
		edge_list = sm.get_edges(l)
		print 'l:', l
		print 'edge_list:', edge_list
		print ''

# test_random_tree(n, repeat_time)


l = [1, [1, 1]]
edge_list = sm.get_edges(l)
print l
print edge_list


###################
# Test tree class #
###################

print 'Test tree class:'

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
t = sm.Tree(edge_list, gp)
# t.print_tree_info()
t.print_node_relation()
print 'root node index:', t.rootNode.index
print 'inorder traversal:', t.in_order_traversal(t.rootNode)
print 'preorder traversal:', t.pre_order_traversal(t.rootNode)

print '\n\n'
print 'test add mutations:'
t.add_mutations_along_edges(t.rootNode)
# t.print_tree_info()

t.print_node_info()









