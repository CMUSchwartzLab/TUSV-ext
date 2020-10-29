
# test cases for combine_copy_nums.py

from combine_copy_nums import *

##############
# Test Cases #
##############
# test case 1:
triplet1_bgns = [0, 1, 4, 7, 8]
triplet1_ends = [0, 3, 6, 7, 17]
triplet1_cps = [3, 2, 3, 2, 1]
triplet1 = [triplet1_bgns, triplet1_ends, triplet1_cps]

triplet2_bgns = [0, 1, 4, 7, 8, 10, 15]
triplet2_ends = [0, 3, 6, 7, 9, 14, 17]
triplet2_cps = [1, 1, 2, 2, 3, 2, 3]
triplet2 = [triplet2_bgns, triplet2_ends, triplet2_cps]

triplet3_bgns = [0, 3, 6]
triplet3_ends = [2, 5, 17]
triplet3_cps = [4, 3, 2]
triplet3 = [triplet3_bgns, triplet3_ends, triplet3_cps]

triplets = [triplet1, triplet2, triplet3]
usages = [0.2, 0.3, 0.5]

[res_bgns, res_ends, res_cps] = combine_copy_nums(triplets, usages)
print 'test case 1:'
print 'res_bgns', res_bgns
print 'res_ends', res_ends
print 'res_cps', res_cps

# tset case 2:
triplet1_bgns = [0,1,4,7,10,14,18]
triplet1_ends = [0,3,6,9,13,17,24]
triplet1_cps = [3, 2, 3, 1, 2, 3, 4]
triplet1 = [triplet1_bgns, triplet1_ends, triplet1_cps]

triplet2_bgns = [0,3,4,10,14,15,18]
triplet2_ends = [2,3,9,13,14,17,24]
triplet2_cps = [1, 1, 2, 2, 2, 3, 1]
triplet2 = [triplet2_bgns, triplet2_ends, triplet2_cps]

triplet3_bgns = [0, 1, 7]
triplet3_ends = [0, 6, 24]
triplet3_cps = [4, 3, 2]
triplet3 = [triplet3_bgns, triplet3_ends, triplet3_cps]

triplets = [triplet1, triplet2, triplet3]
usages = [0.2, 0.7, 0.1]

[res_bgns, res_ends, res_cps] = combine_copy_nums(triplets, usages)
print 'test case 2:'
print 'res_bgns', res_bgns
print 'res_ends', res_ends
print 'res_cps', res_cps

