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


def _calculate_S(S, Pi, U, C, Gam):
    m, l = S.shape
    N, _, _ = Gam.shape
    for p in xrange(0, m):
        for b in xrange(0, l):
            sg_cpnum_est = np.sum([U[p, k] * (Gam[k, b, 0] + Gam[k, b, 1]) for k in xrange(0, N)])
            bp_cpnum_est = np.sum([U[p, k] * C[k, b] for k in xrange(0, N)])
            mod.addConstr(S[p, b] == _get_abs(mod, Pi[p, b] * sg_cpnum_est - bp_cpnum_est))