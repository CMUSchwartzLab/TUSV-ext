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
		while cur != None:  ###xf: go through along the chromosome
			_add_sv_to_dict(svs, others, cur, True, chrm, pm)
			_add_sv_to_dict(svs, others, cur, False, chrm, pm)
			cur = cur.r
		# remove any splits that are not actually breakpoints
		keys_to_remove = []
		for k, v in svs.items():
			if 'mate' not in v:
				keys_to_remove.append(k)
		for k in keys_to_remove:
			del svs[k]
		# add breakpoint copy numbers
		_append_bp_copy_num(svs, others, self.mut)
		return svs, others

	def deepcopy(self):
		c = ChrmProf(self.n, self.chrm, self.pm)
		c.org, muts = _deepcopy_org(self.org)
		muts = sorted(muts, key = lambda x: x.bgn)
		n = len(muts)
		for i in xrange(0, n):
			if i != 0:
				muts[i].l = muts[i-1]
			if i != n-1:
				muts[i].r = muts[i+1]
		c.mut = muts[0]
		return c

	def inv(self, bgn, end):
		if not self._is_in_bounds(bgn, end) or not self._is_splitable(bgn, end):
			return False
		self._2split(bgn, end) # split mutated and original list nodes at bgn and end positions

		self._rev_mut(bgn, end)

		return True

	def rem(self, bgn, end):
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
	def trans(self, from_ChrmProf, ins_Pos, bgn1, end1):
		print(ins_Pos, bgn1, end1)
		if not from_ChrmProf._is_in_bounds(bgn1, end1) or not from_ChrmProf._is_splitable(bgn1, end1):
			return False
		from_ChrmProf._2split(bgn1, end1) # split mutated and original list nodes at bgn and end positions

		insR, head, tail = _copy_from_to(from_ChrmProf.mut, bgn1, end1) ### xf: copied head and tail

		head_, tail_ = _get_head_tail(from_ChrmProf.mut, bgn1, end1) ### xf: original head and tail, to be removed in above codes

		newL_ = head_.l
		newR_ = tail_.r

		if newL_ == None:
			from_ChrmProf.mut = newR_ # change the head of the mut list to right of tail if we are removing head -> tail
		if newL_ != None:
			newL_.r = newR_
		if newR_ != None:
			newR_.l = newL_
		head_.l = None  ### xf: detach the removed segment MutNode
		tail_.r = None

		# remove old nodes from OrgNode children list and delete old nodes
		while head_ != None:
			head_.parent.children.remove(head_) # remove curent MutNode from children list of OrgNode
			prev = head_
			head_ = head_.r
			del prev

		### xf: remove segment finished, start translocation to the new position in current chromosome
		self._split(ins_Pos)
		ins_head = _get_head_ins(self.mut, ins_Pos)
		newL = ins_head.l  ### xf: copy from amp()
		newR = ins_head
		head.l = newL
		tail.r = newR
		if newR != None:
			newR.l = tail
		if newL != None:
			newL.r = head
		# increment bgn and end values for inserted region and segments to right
		seg_diff = ins_Pos - bgn1
		while True:
			head.bgn += seg_diff
			head.end += seg_diff
			if head == tail:
				break
			head = head.r
		head = head.r
		seg_len = end1 - bgn1 + 1
		while head != None:
			head.bgn += seg_len
			head.end += seg_len
			head = head.r
		self.n = self.n + (end1 - bgn1 + 1)
		return from_ChrmProf

	# duplicate region from bgn to end. returns boolean for complete or not
	def amp(self, bgn, end):  	### xf: it uses a linked list data structure, self.mut is always the first MutNode
		if not self._is_in_bounds(bgn, end) or not self._is_splitable(bgn, end):
			return False
		self._2split(bgn, end) # split mutated and original list nodes at bgn and end positions

		insR, head, tail = _copy_from_to(self.mut, bgn, end) # copy list from bgn to end
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
		orgNode2 = orgNode1.split(k)
		for mutNode1 in orgNode1.children:
			mutNode2 = mutNode1.split(k)
			mutNode2.parent = orgNode2
			orgNode2.children.append(mutNode2)

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

	# returns original node that has a mutant at position pos
	def _get_orgNode_mut_pos(pos):
		cur = self.mut
		while cur != None:
			if cur.bgn >= pos and cur.end <= pos:
				return cur.parent
			cur = cur.r
		return None

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


class _OrgNode(_Node): ### xf: OrgNode is a single node which will be splitted during mutation. It will be linked to multiple MutNodes for mapping.
	def __init__(self, bgn, end, chromosome, pm):
		self.children = [] # no mutated sections
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
		r = self.r
		self.r = _MutNode(self.bgn + k, self.end, self.chrm, self.pm, self.is_inv)
		self.r.r = r    # set right of new node to the old node's old right
		self.r.l = self # set left of new node to old node (self)
		if r != None:
			r.l = self.r
		self.end = self.bgn + k - 1
		return self.r

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
def _copy_from_to(head, fm, to):
	oldhead = head
	while head != None and head.bgn != fm: # make head the beginning of where to copy
		head = head.r

	curA = head
	i = 0
	prevB = None
	while curA != None:
		curB = copy.copy(curA)
		curB.parent.children.append(curB) # update parent's children pointers
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

	isLeft = mate.parent.end == matePos
	isAdj = (abs(curPos - matePos) == 1 and curOrgNode.chrm == mateOrgNode.chrm and curOrgNode.pm == mateOrgNode.pm)

	return matePos, isLeft, isAdj, mateOrgNode.chrm, mateOrgNode.pm

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
	if matePos is None:
		return

	curTup = (curPos, curIsLeft, curChrm, curPM)
	mateTup = (matePos, mateIsLeft, mateChrm, matePM)

	
	if curTup not in svs:
		if curChrm == chrm and curPM == pm:
			svs[curTup] = {'total_reads': 0, 'mated_reads': 0}
		else:
			if (curChrm, curPM) in others.keys():
				others[(curChrm, curPM)][curTup] = {'total_reads': 0, 'mated_reads': 0}
			else:
				others[(curChrm, curPM)] = {}
				others[(curChrm, curPM)][curTup] = {'total_reads': 0, 'mated_reads': 0}

	if mateTup not in svs:
		if mateChrm == chrm and matePM == pm:
			svs[mateTup] = {'total_reads': 0, 'mated_reads': 0}
		else:
			if (mateChrm, matePM) in others.keys():
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
				if curTup in others[(curChrm, curPM)] and others[(curChrm, curPM)][curTup]['mate'] == matTup:
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
#
#   DEEP COPY
#

#  input: headA (OrgNode) old head
# output: headB (OrgNode) new head
#         muts (list of MutNode) list of all MutNodes created as children for OrgNodes
def _deepcopy_org(headA):
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
			muts.append(cB)

		curA = curA.r
	return headB, muts
