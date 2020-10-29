import chrm_prof as chpr
import pprint
pp = pprint.PrettyPrinter(indent = 2)

cp = chpr.ChrmProf("AAABBBCCCDDDEEEFFF")
print ''

print cp.amp(4, 10)
cp.pprint()
print ''

print cp.amp(11, 13)
cp.pprint()
print ''

print cp.amp(2, 19)
cp.pprint()
print ''

print cp.rem(0, 6)
cp.pprint()
print ''

print cp.rem(2, 8)
cp.pprint()
print ''

print cp.amp(9, 25)
cp.pprint()
print ''

cp2 = chpr.ChrmProf("RRRRRSSSSSTTTTTUUUUUVVVVV")
print ''

print cp2.amp(5, 14)
cp2.pprint()
print ''

print cp2.inv(3, 20)
cp2.pprint()
print ''

print cp2.inv(15, 30)
cp2.pprint()
print ''

print cp2.rem(10, 16)
cp2.pprint()
print ''

cp3 = chpr.ChrmProf("XXXXYYYYZZZZ")
print ''

print cp3.inv(0, 6)
cp3.pprint()
print ''

print cp3.inv(4, 11)
cp3.pprint()
print ''

# test _get_mated_pos() and _get_cur_pos() functions
cp4 = chpr.ChrmProf("RRRRRSSSSSTTTTTUUUUUVVVVV")
print ''

print cp4.amp(5, 12)
cp4.pprint()
print ''

cur = cp4.mut.r.r
isBgn = True
matePos, isLeft, isAdj = chpr._get_mated_pos(cur, isBgn)

print ''
if isBgn:
	print cur.bgn
else:
	print cur.end
print matePos
print isLeft
print isAdj

print ''
print cp4.inv(3, 25)
cp4.pprint()
print ''

cur = cp4.mut.r.r.r.r
isBgn = True
matePos, isLeft, isAdj = chpr._get_mated_pos(cur, isBgn)

print ''
if isBgn:
	print cur.bgn
else:
	print cur.end
print matePos
print isLeft
print isAdj

curPos, isLeft = chpr._get_cur_pos(cur, isBgn)
print ''
print curPos
print isLeft


# matePos, isLeft, isAdj = chpr._get_mated_pos(cur, False)

# test get_sv_read_nums() function

cov = 20
read_len = 5

cp5 = chpr.ChrmProf("RRRRRSSSSSTTTTTUUUUUVVVVV")
print ''

print cp5.amp(5, 12)
cp5.pprint()
print ''

print cp5.amp(10, 28)
cp5.pprint()
print ''

dic = cp5.get_sv_read_nums(cov, read_len)
pp.pprint(dic)
print ''

# test deepcopy() function

cp5_copied = cp5.deepcopy()
cp5_copied.pprint()
print ''

print cp5_copied.inv(2, 15)
cp5_copied.pprint()
print ''

dic = cp5_copied.get_sv_read_nums(cov, read_len)
pp.pprint(dic)
print ''


dic = cp5.get_sv_read_nums(cov, read_len)
pp.pprint(dic)
print ''