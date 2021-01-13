import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import random
import numpy as np
import multiprocessing as mp

from datetime import datetime
from graphviz import Digraph
from ete2 import Tree          # for creating phylogenetic trees for .xml output
from Bio import Phylo          # for creating phylogenies to export as phylo .xml files
from cStringIO import StringIO # for converting string to file (for creating initial phylo .xml)

# custom modules
sys.path.insert(0, 'model/')
sys.path.insert(0, 'help/')
import solver as sv
import file_manager as fm      # sanitizes file and directory arguments
import generate_matrices as gm # gets F, Q, G, A, H from .vcf files
import printer as pt
import vcf_help as vh

from tusv_1130 import check_valid_input, randomly_remove_segments
from model.solver import np_divide_0


def _calculate_S(Pi, U, C, Gam, m, l):
    S = np.zeros((m, l))
    N, _, _ = Gam.shape
    for p in xrange(0, m):
        for b in xrange(0, l):
            sg_cpnum_est = np.sum([U[p, k] * (Gam[k, b, 0] + Gam[k, b, 1]) for k in xrange(0, N)])
            bp_cpnum_est = np.sum([U[p, k] * C[k, b] for k in xrange(0, N)])
            S[p, b] = np.abs(Pi[p, b] * sg_cpnum_est - bp_cpnum_est)
    return np.sum(S)


def record_true_obj(in_dir, n, c_max, lamb1, lamb2, num_seg_subsamples, should_overide_lambdas):
    F_full, F_phasing_full, Q, G, A, H, bp_attr, cv_attr = gm.get_mats(in_dir)
    print(F_full, Q, G)
    Q, G, A, H, F_full, F_phasing_full = check_valid_input(Q, G, A, H, F_full, F_phasing_full)

    F, F_phasing, Q, org_indxs = randomly_remove_segments(F_full, F_phasing_full, Q, num_seg_subsamples)
    m = len(F)
    l, r = Q.shape
    print(F_phasing.shape)
    # replace lambda1 and lambda2 with input derived values if should_orveride_lamdas was specified
    if should_overide_lambdas:

        lamb1 = float(l + r) / float(l) * float(m) / float(2 * (n - 1))
        lamb2 = float(l + r) / float(l)

    F_seg = F[:, l:].dot(np.transpose(Q))  # [m, l] mixed copy number of segment containing breakpoint
    Pi = np_divide_0(F[:, :l], F_seg)