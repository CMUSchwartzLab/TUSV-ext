import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import random
import numpy as np
import multiprocessing as mp
import pickle

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


def _calculate_S(Pi, U, C, Gam, m, l):
    S = np.zeros((m, l))
    N, _, _ = Gam.shape
    for p in xrange(0, m):
        for b in xrange(0, l):
            sg_cpnum_est = np.sum([U[p, k] * (Gam[k, b, 0] + Gam[k, b, 1]) for k in xrange(0, N)])
            bp_cpnum_est = np.sum([U[p, k] * C[k, b] for k in xrange(0, N)])
            S[p, b] = np.abs(Pi[p, b] * sg_cpnum_est - bp_cpnum_est)
    return S


def _calculate_Gamma(Q, C, n):
    l, r = Q.shape
    N = 2 * n - 1
    Gam = np.zeros((N, l, 2))
    for k in xrange(0, N):
        for b in xrange(0, l):  # define copy num of segment containing breakpoint
            Gam[k, b, 0] = np.sum([Q[b, s] * C[k, l + s] for s in xrange(0, r)]) ### xf: change to new Gamma and C matrix
            Gam[k, b, 1] = np.sum([Q[b, s] * C[k, l + s + r] for s in xrange(0, r)])
    return Gam


def _calculate_R(C, edge_list, l):
    N = C.shape[0]
    R = np.zeros((N, N))
    for edge in edge_list:
        p = edge[0]-1
        c = edge[1]-1
        R[p, c] = np.sum(np.abs(C[p,l:] - C[c,l:]))
    return R


def _calculate_obj_val(F_phasing, C, U, R, S, lambda1, lambda2):
    F_hat = np.matmul(U, C)
    obj_val = np.sum(np.abs(F_hat - F_phasing)) + lambda1 * np.sum(R) + lambda2*np.sum(S)
    return obj_val


def np_divide_0(a, b):
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        c[~ np.isfinite(c)] = 0  # -inf inf NaN
    return c