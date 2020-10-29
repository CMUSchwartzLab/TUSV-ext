
# test cases for combine_copy_nums.py

from combine_copy_nums import *

##############
# Test Cases #
##############
# test case 1:
triplet1_bgns = [1, 31, 61]
triplet1_ends = [30, 60, 100]
triplet1_cps = [3, 2, 3]
triplet1 = [triplet1_bgns, triplet1_ends, triplet1_cps]

triplet2_bgns = [1, 21, 51]
triplet2_ends = [20, 50, 100]
triplet2_cps = [1, 1, 2]
triplet2 = [triplet2_bgns, triplet2_ends, triplet2_cps]

triplet3_bgns = [1, 41, 81]
triplet3_ends = [40, 80, 100]
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
triplet1_bgns = [1, 31, 61]
triplet1_ends = [30, 60, 100]
triplet1_cps = [3, 2, 3]
triplet1 = [triplet1_bgns, triplet1_ends, triplet1_cps]

triplet2_bgns = [1, 21, 60]
triplet2_ends = [20, 59, 100]
triplet2_cps = [1, 1, 2]
triplet2 = [triplet2_bgns, triplet2_ends, triplet2_cps]

triplets = [triplet1, triplet2]
usages = [0.2, 0.8]

[res_bgns, res_ends, res_cps] = combine_copy_nums(triplets, usages)
print 'test case 2:'
print 'res_bgns', res_bgns
print 'res_ends', res_ends
print 'res_cps', res_cps

