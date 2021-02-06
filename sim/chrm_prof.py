#     file: chrm_prof.py
#   author: Jesse Eaton
#  created: 9/27/2017
# modified: 9/30/2017
#  purpose: ChrmProf class. ChrmProf is a chromosome profile and contains a mutated chrm with
#             references to the mutations' original position. Keeps track of copy numbers.

# imports
import sys
import copy
import math
import random
import numpy as np

# helpers
def printnow(s):
	s = str(s)
	sys.stdout.write(s)
	sys.stdout.flush()
def printerr(s):
	sys.stderr.write(s)
	sys.stderr.flush()

class ChrmProf:  ### xf: the profile specifically for one chromosome (allele specific)
	# n (length of chromosome)
	def __init__(self, n, chromosome, pm):
		self.n = n
		self.chrm = chromosome
		self.pm = pm
		self.org = _OrgNode(0, n - 1, chromosome, pm)
		self.mut = _MutNode(0, n - 1, chromosome, pm)
		self.org.children.append(self.mut)
		self.mut.parent = self.org

	# output: bgns (list of int) [n] beginning positions for each segment
	#         ends (list of int) [n] ending positions for each segment
	#         cps  (list of int) [n] copy number for each segment
	def get_copy_nums(self):
		cur = self.org
		bgns, ends, cps = [], [], []
		while cur != None:
			bgns.append(cur.bgn)
			ends.append(cur.end)
			cps.append(len(cur.children))
			cur = cur.r
		return bgns, ends, cps

	def get_sv_read_nums(self, cov, read_len, chrm, pm):
		n = self.n
		svs = {}
		others = {}
		cur = self.mut
		#print(self.chrm,self.pm)
		while cur != None:  ###xf: go through along the chromosome
			svs, others = _add_sv_to_dict(svs, others, cur, True, chrm, pm)
			svs, others = _add_sv_to_dict(svs, others, cur, False, chrm, pm)
			#print(cur.bgn,cur.end, _get_org_pos(cur,True)[0], _get_org_pos(cur,False)[0])
			cur = cur.r

		# remove any splits that are not actually breakpoints
		keys_to_remove = []
		for k, v in svs.items():
			if 'mate' not in v:
				keys_to_remove.append(k)
		for k in keys_to_remove:
			del svs[k]
		for k, v in others.items():
			keys_to_remove = []
			for k1, v1 in others[k].items():
				if 'mate' not in v1:
					keys_to_remove.append(k1)
			for k2 in keys_to_remove:
				del others[k][k2]

		# add breakpoint copy numbers
		svs, others = _append_bp_copy_num(svs, others, self.mut)
		return svs, others

	def get_snvs(self, snvs):
		cur = self.mut
		while cur != None:  ###xf: go through along the chromosome
			for snv_mut in cur.SNV_Mut_children:
				snv_org = snv_mut.SNV_Org_parent
				snv_tuple = (snv_org.chrm, snv_org.pos)
				if snv_tuple not in snvs.keys():
					snvs[snv_tuple] = {'cn': 0}
				snvs[snv_tuple]['cn'] += 1
			cur = cur.r
		return snvs

	### xf: add SNVs
	def point_mutation(self, pos):
		splitMut = self.mut
		while splitMut != None and not (splitMut.bgn <= pos and pos <= splitMut.end):
			splitMut = splitMut.r

		SNV_MutNode = _SNV_MutNode(pos, self.chrm, self.pm)

		SNV_MutNode.Mut_parent = splitMut
		splitMut.SNV_Mut_children.append(SNV_MutNode)

		orgNode = splitMut.parent
		if splitMut.is_inv:
			org_pos = orgNode.end - (pos - splitMut.bgn)
		else:
			org_pos = pos
		SNV_OrgNode = _SNV_OrgNode(org_pos, orgNode.chrm, orgNode.pm)

		SNV_OrgNode.Org_parent = orgNode
		orgNode.SNV_Org_children.append(SNV_OrgNode)

		SNV_OrgNode.SNV_Mut_children.append(SNV_MutNode)
		SNV_MutNode.SNV_Org_parent = SNV_OrgNode
		return True

	def deepcopy_(self, other_muts):
		c = ChrmProf(self.n, self.chrm, self.pm)
		c.org, muts, other_muts = _deepcopy_org(self.org, other_muts)
		return c, muts, other_muts

	def inv(self, bgn, end, snv):
		print("inv",self.chrm, self.pm, bgn,end)
		if not self._is_in_bounds(bgn, end) or not self._is_splitable(bgn, end):
			return False
		if bgn == end:
			return False
		self._2split(bgn, end) # split mutated and original list nodes at bgn and end positions

		self._rev_mut(bgn, end)
		if snv:
			self._rev_mut_snv(bgn, end)

		return True

	def rem(self, bgn, end, snv):
		print("rem",self.chrm, self.pm, bgn, end)
		if not self._is_in_bounds(bgn, end) or not self._is_splitable(bgn, end):
			return False
		self._2split(bgn, end) # split mutated and original list nodes at bgn and end positions

		head, tail = _get_head_tail(self.mut, bgn, end) ### xf: "# head is MutNode with bgn and tail is MutNode with end."

		newL = head.l
		newR = tail.r

		if newL == None:
			self.mut = newR # change the head of the mut list to right of tail if we are removing head -> tail
		if newL != None:
			newL.r = newR
		if newR != None:
			newR.l = newL
		head.l = None  ### xf: detach the removed segment MutNode
		tail.r = None

		# remove old nodes from OrgNode children list and delete old nodes
		while head != None:
			head.parent.children.remove(head) # remove curent MutNode from children list of OrgNode
			if snv:
				for snv_child in head.SNV_Mut_children:
					snv_child.SNV_Org_parent.remove(snv_child)
					snv_child.Mut_parent.remove(snv_child)
					del snv_child
			prev = head
			head = head.r
			del prev

		# decrement bgn and end values for segments to right of deleted region
		seg_len = end - bgn + 1
		cur = newR
		while cur != None:
			cur.bgn -= seg_len
			cur.end -= seg_len
			cur = cur.r

		self.n = self.n - (end - bgn + 1)

		return True

	###xf: add translocation
	def trans(self, from_ChrmProf, ins_Pos, bgn1, end1, snv):
		print('trans', self.chrm + ' ' + str(self.pm), ins_Pos, from_ChrmProf.chrm + ' ' + str(from_ChrmProf.pm), bgn1, end1)
		if not from_ChrmProf._is_in_bounds(bgn1, end1) or not from_ChrmProf._is_splitable(bgn1, end1):
			return False
		if not self._is_in_bounds(ins_Pos, ins_Pos) or not self._is_splitable_one(ins_Pos):
			return False
		from_ChrmProf._2split(bgn1, end1) # split mutated and original list nodes at bgn and end positions

		#insR, head, tail = _copy_from_to(from_ChrmProf.mut, bgn1, end1) ### xf: copied head and tail
		#print(_get_org_pos(head, True)[0])
		head_, tail_ = _get_head_tail(from_ChrmProf.mut, bgn1, end1) ### xf: original head and tail, to be removed in above codes

		newL_ = head_.l
		newR_ = tail_.r

		if newL_ == None:
			from_ChrmProf.mut = newL_ # change the head of the mut list to right of tail if we are removing head -> tail
		if newL_ != None:
			newL_.r = newR_
		if newR_ != None:
			newR_.l = newL_
		# head_.l = None  ### xf: detach the removed segment MutNode
		# tail_.r = None

		# remove old nodes from OrgNode children list and delete old nodes
		# while head_ != None:
		# 	head_.parent.children.remove(head_) # remove curent MutNode from children list of OrgNode
		# 	prev = head_
		# 	head_ = head_.r
		# 	del prev

		seg_len = end1 - bgn1 + 1
		right = tail_.r


		### xf: remove segment finished, start translocation to the new position in current chromosome
		self._split(ins_Pos)
		ins_head = _get_head_ins(self.mut, ins_Pos)
		newL = ins_head.l  ### xf: copy from amp()
		newR = ins_head
		head_.l = newL
		tail_.r = newR
		if newR != None:
			newR.l = tail_
		if newL != None:
			newL.r = head_
		# increment bgn and end values for inserted region and segments to right
		# decrement bgn and end values for translocated segment end to right
		while right is not None:
			right.bgn -= seg_len
			right.end -= seg_len
			if snv:
				for snv_child in right.SNV_Mut_children:
					snv_child.pos -= seg_len
			right = right.r
		# increment bgn and end values for inserted region and segments to right
		seg_diff = ins_Pos - bgn1
		head_arc = head_
		while True:
			head_.bgn += seg_diff
			head_.end += seg_diff
			head_.chrm = self.chrm
			head_.pm = self.pm
			if snv:
				for snv_child in head_.SNV_Mut_children:
					snv_child.pos += seg_diff
					snv_child.chrm = self.chrm
					snv_child.pm = self.pm
			if head_ == tail_:
				break
			head_ = head_.r
		head_ = head_.r


		while head_ != None:
			head_.bgn += seg_len
			head_.end += seg_len
			if snv:
				for snv_child in head_.SNV_Mut_children:
					snv_child.pos += seg_len
			head_ = head_.r
		self.n = self.n + (end1 - bgn1 + 1)
		from_ChrmProf.n -= (end1 - bgn1 + 1)
		# mut = from_ChrmProf.mut
		# while mut is not None:
		# 	print(mut.bgn, mut.end)
		# 	mut = mut.r
		# mut2 = self.mut
		# while mut2 is not None:
		# 	print(mut2.bgn, mut2.end)
		# 	mut2 = mut2.r
		# print(from_ChrmProf.mut.bgn, self.mut.bgn)
		return from_ChrmProf

	# duplicate region from bgn to end. returns boolean for complete or not
	def amp(self, bgn, end, snv):  	### xf: it uses a linked list data structure, self.mut is always the first MutNode
		print("amp",self.chrm, self.pm, bgn, end)
		if not self._is_in_bounds(bgn, end) or not self._is_splitable(bgn, end):
			return False
		self._2split(bgn, end) # split mutated and original list nodes at bgn and end positions

		insR, head, tail = _copy_from_to(self.mut, bgn, end, snv) # copy list from bgn to end
		### xf: copy means copy the whole identity including the parent-children relationship
		### xf: duplicate two consecutive nodes, to visualize it looks like: insR-head-.....-tail-insL for MutNode
		insL = insR.r # node to go after tail
		insR.r = head
		head.l = insR
		tail.r = insL
		if insL != None:
			insL.l = tail

		# increment bgn and end values for inserted region and segments to right
		seg_len = end - bgn + 1
		while head != None:
			head.bgn += seg_len
			head.end += seg_len
			if snv:
				for snv_child in head.SNV_Mut_children:
					snv_child.pos += seg_len
			head = head.r

		self.n = self.n + (end - bgn + 1)
		return True

	# split bgn and end positions if needed. do not need to split at start or terminal of chromosome
	def _2split(self, bgn, end):
		n = self.n
		if bgn > 0:
			self._split(bgn)     # do not split if bgn is 0 (start of chrm)
		if end + 1 < n:
			self._split(end + 1) # do not split if end is n-1 (end of chrm)

	# splits node at position k in mut and org. k is bgn of right offspring
	def _split(self, k):

		# return if k out of bounds
		if k == 0 or k >= self.n:
			printerr('cannot split at k = ' + str(k) + '\n')
			return

		# find orgNode along the genome corresponding to the mutNode where split will occur
		splitMut = self.mut
		while splitMut != None and not (splitMut.bgn <= k and k <= splitMut.end):
			splitMut = splitMut.r
		orgNode1 = splitMut.parent

		if splitMut.bgn == k or splitMut.end == k-1: # should not split b/c this was already split
			return

		k = k - splitMut.bgn # make k the new length of the old node (node that will be split)
		                     # is is now the number of nucletides in from a node where it should be split

		# split orgNode1 and all its children
		### xf: split the orgNode with different position according to whether this mutNode is inv or not
		if splitMut.is_inv:
			orgNode2 = orgNode1.split(splitMut.end - k - splitMut.bgn + 1)
			### xf: relink the parent-children relationship for SVs and SNVs
			temp_children_list = orgNode1.SNV_Org_children.copy()
			for child in temp_children_list:
				if child.pos >= splitMut.end - k - splitMut.bgn + 1:
					child.Org_parent = orgNode2
					orgNode2.SNV_Org_children.append(child)
					orgNode1.SNV_Org_children.remove(child)
		else:
			orgNode2 = orgNode1.split(k)
			temp_children_list = orgNode1.SNV_Org_children.copy()
			for child in temp_children_list:
				if child.pos >= k:
					child.Org_parent = orgNode2
					orgNode2.SNV_Org_children.append(child)
					orgNode1.SNV_Org_children.remove(child)
		for mutNode1 in orgNode1.children:
			mutNode2 = mutNode1.split(k)
			mutNode2.parent = orgNode2
			orgNode2.children.append(mutNode2)
			temp_children_list = mutNode1.SNV_Mut_children
			if splitMut.is_inv:
				for mutchild in temp_children_list:
					if mutchild.pos <= splitMut.bgn + k - 1:
						mutchild.Mut_parent = mutNode2
						mutNode2.SNV_mut_children.append(mutchild)
						mutNode1.SNV_Mut_children.remove(mutchild)
			else:
				for mutchild in temp_children_list:
					if mutchild.pos >= splitMut.bgn + k:
						mutchild.Mut_parent = mutNode2
						mutNode2.SNV_Mut_children.append(mutchild)
						mutNode1.SNV_Mut_children.remove(mutchild)

	def _is_in_bounds(self, bgn, end):
		n = self.n
		if bgn < 0 or bgn > end or end >= n:
			printerr('chromosome has length ' + str(n) + '. cannot mutate [' + str(bgn) + ', ' + str(end) + ']')
			return False
		return True

	# returns True if bgn and end do not match any positions already in mutated list
	def _is_splitable(self, bgn, end):
		n = self.n
		if bgn != 0 and _is_already_mut_bgn(self.mut, bgn):
			return False
		if end + 1 < n and _is_already_mut_end(self.mut, end):
			return False
		return True

	def _is_splitable_one(self, insPos):
		n = self.n
		if insPos != 0 and _is_already_mut_bgn(self.mut, insPos):
			return False
		return True

	# returns original node that has a mutant at position pos
	def _get_orgNode_mut_pos(pos):
		cur = self.mut
		while cur != None:
			if cur.bgn >= pos and cur.end <= pos:
				return cur.parent
			cur = cur.r
		return None

	def _rev_mut_snv(self, bgn, end):
		cur = self.mut
		while cur != None and cur.bgn != bgn: # node with bgn should exist
			cur = cur.r
		for snvMut in cur.SNV_Mut_children:
			snvMut.pos = bgn + end - snvMut.pos

	# reverses doubly linked list starting from node with bgn of bgn to node with end of end
	def _rev_mut(self, bgn, end):
		cur = self.mut
		while cur != None and cur.bgn != bgn: # node with bgn should exist
			cur = cur.r
		ih = cur   # inner head
		oh = cur.l # outer head
		while cur != None and cur.end != end: # node with end should exist
			cur = cur.r
		it = cur   # inner tail
		ot = cur.r # outer tail

		# set region bgn and end for calculating new positions of segments
		rgbgn, rgend = ih.bgn, it.end
		
		cur = ih
		while cur != ot: # reverse linked list ih -> ... -> it
			prv = cur.l
			nxt = cur.r
			cur.l = nxt
			cur.r = prv
			pos_diff = _get_inv_pos_diff(rgbgn, rgend, cur.bgn, cur.end)
			cur.bgn += pos_diff
			cur.end += pos_diff
			cur.is_inv = not cur.is_inv
			cur = nxt

		it.l = oh      # connect head and tail of internally reversed linked list to outer list
		ih.r = ot
		if oh != None:
			oh.r = it
		if ot != None:
			ot.l = ih

		if oh == None:
			self.mut = it # set head of mut if we inverted start of list

	def pprint(self):
		printnow('lists:\n')
		for lst_name, cur in {'org': self.org, 'mut': self.mut}.iteritems():
			printnow(lst_name + ': ')
			while cur != None:
				cur.pprint()
				cur = cur.r
			printnow('\n')
		printnow('relations:\n')
		cur = self.org
		while cur != None:
			kid_pos_strs = [ kid.get_pos_str() for kid in cur.children ]
			printnow('[' + str(cur.bgn) + ',' + str(cur.end) + '] -> ' + ', '.join(kid_pos_strs) + '\n')
			cur = cur.r
		printnow('copy numbers: ' + str(self.get_copy_nums()[2]) + '\n')

class _Node:
	def pprint(self):
		s = self.get_pos_str()
		if self.r != None:
			s += '->'
		else:
			s += '-v'
		printnow(s)
	def get_pos_str(self):
		return '[' + str(self.bgn) + ',' + str(self.end) + ']'

### xf: add class for SNV nodes:
class _SNV_MutNode:
	def __init__(self, pos, chromosome, pm):
		self.Mut_parent = None
		self.SNV_Org_parent = None
		self.pos = pos
		self.chrm = chromosome
		self.pm = pm

	def copy(self):
		return _SNV_MutNode(self.pos, self.chrm, self.pos)

class _SNV_OrgNode:
	def __init__(self, pos, chromosome, pm):
		self.Org_parent = None
		self.SNV_Mut_children = []
		self.pos = pos
		self.chrm = chromosome
		self.pm = pm

	def copy(self):
		return _SNV_OrgNode(self.pos, self.chrm, self.pos)


class _OrgNode(_Node): ### xf: OrgNode is a single node which will be splitted during mutation. It will be linked to multiple MutNodes for mapping.
	def __init__(self, bgn, end, chromosome, pm):
		self.children = [] # no mutated sections
		self.SNV_Org_children = []
		self.l = None      # no left or right pointers
		self.r = None
		self.bgn = bgn
		self.end = end
		self.chrm = chromosome
		self.pm = pm

	# returns pointer to new sibling on right. k (int) means k + self.begin is bgn of new sibling
	def split(self, k):
		r = self.r
		self.r = _OrgNode(self.bgn + k, self.end, self.chrm, self.pm)
		self.r.r = r    # set right of new node to the old node's old right
		self.r.l = self # set left of new node to old node (self)
		if r != None:
			r.l = self.r
		self.end = self.bgn + k - 1
		return self.r

	def copy(self):
		return _OrgNode(self.bgn, self.end, self.chrm, self.pm)


class _MutNode(_Node):
	def __init__(self, bgn, end, chromosome, pm, is_inv = False):
		self.parent = None
		self.SNV_Mut_children = []
		self.l = None
		self.r = None
		self.bgn = bgn
		self.end = end
		self.is_inv = is_inv
		self.chrm = chromosome
		self.pm = pm

	def copy(self):
		return _MutNode(self.bgn, self.end, self.chrm, self.pm, self.is_inv)

	# returns pointer to new sibling on right. k (int) means k + self.begin is bgn of new sibling
	def split(self, k):
		if not self.is_inv:
			r = self.r
			self.r = _MutNode(self.bgn + k, self.end, self.chrm, self.pm, self.is_inv)
			self.r.r = r    # set right of new node to the old node's old right
			self.r.l = self # set left of new node to old node (self)
			if r != None:
				r.l = self.r
			self.end = self.bgn + k - 1
			return self.r
		### xf: add the scenario where the split point might be within an inversed segment
		else:
			l = self.l
			self.l = _MutNode(self.bgn, self.bgn + k - 1, self.chrm, self.pm, self.is_inv)
			self.l.l = l
			self.l.r = self
			if l != None:
				l.r = self.l
			self.bgn = self.bgn + k
			return self.l

	def get_pos_str(self):
		if self.is_inv:
			return 'i[' + str(self.bgn) + ',' + str(self.end) + ']'
		return '[' + str(self.bgn) + ',' + str(self.end) + ']'

# helpers

# head is MutNode with bgn and tail is MutNode with end. must _2split() before so these exist!
def _get_head_tail(cur, bgn, end):
	while cur != None and cur.bgn != bgn:
		cur = cur.r
	head = cur
	while cur != None and cur.end != end:
		cur = cur.r
	tail = cur
	return head, tail

def _get_head_ins(cur, insPos):
	while cur != None and cur.bgn != insPos:
		cur = cur.r
	head = cur
	return head

# cur (mutNode). returns True if bgn is already in mutNode list
def _is_already_mut_bgn(cur, bgn):
	while cur != None:
		if bgn == cur.bgn:
			return True
		cur = cur.r
	return False
def _is_already_mut_end(cur, end):
	while cur != None:
		if end == cur.end:
			return True
		cur = cur.r
	return False

# fm (int) is bgn index of one of the nodes. to (int) is end of one of the nodes
### xf: head means the current start MutNode
def _copy_from_to(head, fm, to, snv):
	oldhead = head
	while head != None and head.bgn != fm: # make head the beginning of where to copy
		head = head.r

	curA = head
	i = 0
	prevB = None
	while curA != None:
		curB = copy.copy(curA)
		curB.parent.children.append(curB) # update parent's children pointers
		if snv:
			for snv_child in curA.SNV_Mut_children:
				mutB = copy.copy(snv_child)
				mutB.SNV_Org_parent = snv_child.SNV_Org_parent
				snv_child.SNV_Org_parent.SNV_Mut_children.append(mutB)
				curB.SNV_Mut_children.append(mutB)
				mutB.Mut_parent = snv_child.Mut_parent
		if i == 0:
			headB = curB
			i += 1
		curB.l = prevB
		if prevB != None:
			prevB.r = curB
		curB.r = None
		prevB = curB
		if curA.end == to:
			return curA, headB, curB
		curA = curA.r
	printerr('should not get here')
	cur = oldhead
	printnow('\n\n')
	while cur != None:
		cur.pprint()
		cur = cur.r
	printnow('\n\n')

#  input: rgbgn (int) region beginning. position of bgn for first node of list to be inverted
#         rgend (int) region ending.    position of end for last  node of list to be inverted
#         sgbgn (int) segment beginning. bgn of current segment being flipped
#         sgend (int) segment ending.    end of current segment being flipped
# output: d (int) change in position to be applied to bgn and end of current node
def _get_inv_pos_diff(rgbgn, rgend, sgbgn, sgend):
	rglen = rgend - rgbgn + 1                            # region length
	sglen = sgend - sgbgn + 1                            # segment length
	mid = int(math.ceil(float(rglen) / 2.0) + rgbgn - 1) # position of midpoint of region
	l = mid - sgend                                      # distance from midpoint to segment ending
	if rglen % 2 == 0: # even
		return 2 * l + sglen
	return 2 * l + sglen - 1

def tri_split_str(s, bgn, end):
	s1 = s[:bgn]
	s2 = s[bgn:end+1]
	s3 = s[end+1:]
	return s1, s2, s3

# returns True if x is between a and b inclusive
def _is_between(x, a, b):
	return (a <= x) and (x <= b)

#  input: cur (MutNode) current mutant node the read maps to
#         isBgn (bool) True if read mapped to cur.bgn. False if mapped to cur.end
# output: pos (int) mate OrgNode's bgn or end position (depends on isBgn). this is position of mate bp
#         isLeft (bool) mate OrgNode's orientation. True if cur was found next to end on mate. False if
#                                                           cur was found next to bgn on mate
#         isAdj (bool) True if mate is adjacent to cur in original genome. this will not increment num mated reads
def _get_mated_pos(cur, isBgn): ### xf: find out any breakpoint 123|124
	mate = cur.r  ### xf: mated pairs are doomed to be adjacent in MutNode,  mate | cur or cur | mate
	if isBgn:
		mate = cur.l

	if mate == None:
		return None, None, None, None, None

	curPos, curOrgNode = _get_org_pos(cur, isBgn)
	matePos, mateOrgNode = _get_org_pos(mate, not isBgn)

	isLeft = isBgn
	if cur.is_inv:
		isLeft = not isBgn
	isAdj = (abs(curPos - matePos) == 1 and curOrgNode.chrm == mateOrgNode.chrm and curOrgNode.pm == mateOrgNode.pm)

	return matePos, isLeft, isAdj, mateOrgNode.chrm, mateOrgNode.pm


def _get_org_pos_snv(snv_MutNode):
	return snv_MutNode.SNV_Org_parent.pos


# node (MuteNode), isBgn (bool) True if considering left pos on mutant. returns position of org node
def _get_org_pos(node, isBgn):
	if node.is_inv:
		isBgn = not isBgn
	if isBgn:
		return node.parent.bgn, node.parent
	return node.parent.end, node.parent

def _get_cur_pos(cur, isBgn):
	oCur = cur.parent
	isLeft = not isBgn
	if cur.is_inv:
		isLeft = not isLeft
	if isLeft:
		return oCur.end, isLeft, oCur.chrm, oCur.pm
	return oCur.bgn, isLeft, oCur.chrm, oCur.pm


def _add_sv_to_dict(svs, others, cur, isBgn, chrm, pm):
	matePos, mateIsLeft, isAdj, mateChrm, matePM = _get_mated_pos(cur, isBgn)   ### xf: mate | cur
	curPos, curIsLeft, curChrm, curPM = _get_cur_pos(cur, isBgn)
	curTup = (curPos, curIsLeft, curChrm, curPM)
	mateTup = (matePos, mateIsLeft, mateChrm, matePM)
	#print(curTup, mateTup)
	if matePos is None:
		return svs, others

	# curTup = (curPos, curIsLeft, curChrm, curPM)
	# mateTup = (matePos, mateIsLeft, mateChrm, matePM)

	if curTup not in svs:
		if curChrm == chrm and curPM == pm:
			svs[curTup] = {'total_reads': 0, 'mated_reads': 0}
		else:
			if (curChrm, curPM) in others.keys():
				if curTup not in others[(curChrm, curPM)].keys():
					others[(curChrm, curPM)][curTup] = {'total_reads': 0, 'mated_reads': 0}
			else:
				others[(curChrm, curPM)] = {}
				others[(curChrm, curPM)][curTup] = {'total_reads': 0, 'mated_reads': 0}

	if mateTup not in svs:
		if mateChrm == chrm and matePM == pm:
			svs[mateTup] = {'total_reads': 0, 'mated_reads': 0}
		else:
			if (mateChrm, matePM) in others.keys():
				if mateTup not in others[(mateChrm, matePM)].keys():
					others[(mateChrm, matePM)][mateTup] = {'total_reads': 0, 'mated_reads': 0}
			else:
				others[(mateChrm, matePM)] = {}
				others[(mateChrm, matePM)][mateTup] = {'total_reads': 0, 'mated_reads': 0}

	if curChrm == chrm and curPM == pm:
		svs[curTup]['total_reads'] += 1
		if not isAdj:  ### xf: not include normal adjacent pairs (need to adjust if we include SNVs)
			svs[curTup]['mated_reads'] += 1
			svs[curTup]['mate'] = mateTup
			if mateChrm == chrm and matePM == pm:
				svs[mateTup]['mate'] = curTup
			else:
				others[(mateChrm, matePM)][mateTup]['mate'] = curTup
	else:
		others[(curChrm, curPM)][curTup]['total_reads'] += 1
		if not isAdj: ### xf: not include normal adjacent pairs (need to adjust if we include SNVs)
			others[(curChrm, curPM)][curTup]['mated_reads'] += 1
			others[(curChrm, curPM)][curTup]['mate'] = mateTup
			if mateChrm == chrm and matePM == pm:
				svs[mateTup]['mate'] = curTup
			else:
				others[(mateChrm, matePM)][mateTup]['mate'] = curTup
	return svs, others



def _append_bp_copy_num(svs, others, mut_head):
	cur = mut_head
	while cur != None:
		for isBgn in [True, False]:
			curPos, curIsLeft, curChrm, curPM = _get_cur_pos(cur, isBgn)
			matPos, matIsLeft, _, mateChrm, matePM = _get_mated_pos(cur, isBgn)
			curTup = (curPos, curIsLeft, curChrm, curPM)
			matTup = (matPos, matIsLeft, mateChrm, matePM)
			if curTup in svs and svs[curTup]['mate'] == matTup:
				if 'copy_num' not in svs[curTup]:
					svs[curTup]['copy_num'] = 0
				svs[curTup]['copy_num'] += 1
			elif (curChrm, curPM) in others.keys():
				if curTup in others[(curChrm, curPM)]:
					if others[(curChrm, curPM)][curTup]['mate'] == matTup:
						if 'copy_num' not in others[(curChrm, curPM)][curTup]:
							others[(curChrm, curPM)][curTup]['copy_num'] = 0
						others[(curChrm, curPM)][curTup]['copy_num'] += 1
		cur = cur.r

	# add copy number of zeros for breakpoints that were deleted
	for tup, val in svs.iteritems():
		if 'copy_num' not in val:
			svs[tup]['copy_num'] = 0
	for tup1, others_chr in others.iteritems():
		for tup2, val in others_chr.iteritems():
			if 'copy_num' not in val:
				others[tup1][tup2]['copy_num'] = 0
	#print("svs after", svs, "others after", others)
	return svs, others
#   DEEP COPY
#

#  input: headA (OrgNode) old head
# output: headB (OrgNode) new head
#         muts (list of MutNode) list of all MutNodes created as children for OrgNodes
def _deepcopy_org(headA, other_muts):
	chrm = headA.chrm
	pm = headA.pm
	curA = headA
	i = 0
	prvB = None
	muts = []
	while curA != None:
		curB = curA.copy()
		if i == 0:
			headB = curB
			i += 1
		curB.l = prvB
		if prvB != None:
			prvB.r = curB
		prvB = curB

		# create all mut children
		for cA in curA.children:
			cB = cA.copy()
			cB.parent = curB
			curB.children.append(cB)
			if cA.chrm == chrm and cA.pm == pm:
				muts.append(cB)
			else:
				if (cA.chrm, cA.pm) not in other_muts.keys():
					other_muts[(cA.chrm, cA.pm)] = []
				other_muts[(cA.chrm, cA.pm)].append(cB)
		curA = curA.r
	return headB, muts, other_muts

###xf: deepcopy mutNode
def _deepcopy_mut(headA):
	curA = headA
	i = 0
	prvB = None
	muts = []
	while curA != None:
		curB = curA.copy()
		if i == 0:
			headB = curB
			i += 1
		curB.l = prvB
		if prvB != None:
			prvB.r = curB
		prvB = curB

		# create all mut children
		pA = curA.parent
		pB = pA.copy()


		curA = curA.r
	return headB, muts