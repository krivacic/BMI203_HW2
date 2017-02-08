import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically, compute_similarity
import glob
import os
from .utils import Atom, Residue, ActiveSite
import glob
import os
import sys
import random
import numpy as np
from shutil import copyfile
from .io import read_active_sites, read_active_site, write_clustering, write_mult_clusterings


# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 5:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file> <number of iterations>")
    sys.exit(0)

active_sites = read_active_sites(sys.argv[2])


# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(active_sites)
    write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clusterings = cluster_hierarchically(active_sites)
    write_mult_clusterings(sys.argv[3], clusterings)

#compute_similarity()
if sys.argv[1][0:2] == '-S':
    print("Testing similarity metric")
    """ This script takes every file in 'data' and runs it through compute_similarity with
    every other file, without being redundant (i.e. each pair only goes through
    compute_similarity once) """

    #for i in active_sites:
        #matrix_create(i)
    score_list = []
    for i in active_sites:
        site1 = i
        for j in active_sites:
            if j not in score_list:
                print(i,j)
                site2 = j
                score1 = compute_similarity(site1,site2)
                print(score1)
                #print("site1",site1,"site2", site2)
        score_list.append(i)

#site1 = read_active_site("data/13052.pdb")
#site2 = read_active_site("data/81816.pdb")
#score = compute_similarity(site1,site1)
#print(site1.residues[1].atoms[1].coords)
#print(site2.residues[1].atoms[1].coords)
#print(score)
