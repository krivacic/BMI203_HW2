from .utils import Atom, Residue, ActiveSite
import glob
import os
import sys
import random
import numpy as np
from shutil import copyfile
from .io import read_active_sites, read_active_site, write_clustering, write_mult_clusterings

#active_sites = read_active_sites(sys.argv[2])
#active_sites = read_active_sites(sys.argv[2])
"""
Usage:
From parent folder:
python -m hw2skeleton [-P | -H | -S] <path to data> <path to output file.txt> #iterations #clusters(if -P)


Functions:

silhouette_P(clusters): computes silhouette score for the partitioning algorithm.

silhouette_H(clusters): computes silhouette score for the hierarchical algorithm.

centroid(site_a): translates active site so its centroid is 0,0,0.

matrix_create(site_a,site_b): creates a matrix of atom coordinates from active
site a and active site b. Since the rotation matrix only works when the dimensions
of these matrixes are the same, the function also calls an alignment function
in order to determine which residues/atoms to trim from the longer of the
two active sites.

find_rotation_matrix(P,Q): finds a matrix U that optimally rotates matrix P
onto matrix Q, where matrices P and Q are matrices of atom coordinates.

rotate_matrix(P,Q): rotates matrix P onto matrix Q using the rotation matrix U.

res_similarity(aa1,aa2): quick function to determine the similarity between two
residues, using the BLOSUM62 amino acid substitution matrix.

sequence_alignment(residues_a,residues_b): performs a sequence alignment of residue
lists from two active sites, which helps determine which residues to leave out in
the calculation of RMSD.

compute_similarity(site_a,site_b): computes the RMSD of all rotated atoms. This
function itself does not take into account residue similarity, but the impact of
this should be minimal since the sequence alignment already accounts for residue
similarity. Essentially, the similarity score is the RMSD of the most similar residues
in the active site.

find_distances(active_sites): creates a dictionary with all distances so that they only
need to be computed once.

clustering_by_partitioning(): clusters all active sites via K-medioid partitioning
method.

cluster_hierarchically(): cluster all active sites hierarically. Uses agglomorative
clustering with linkage by averages.


"""

def silhouette_P(clusters,distances):

    dist = distances
    Clist = []
    for key in clusters:
        sitelist = []
        for subkey in clusters[key]:
            site1 = subkey
            sitelist.append(subkey)
            ilist = []
            for subkey2 in clusters[key]:
                if subkey2 != site1:
                    site2 = subkey2
                    D = dist[site1][site2]
                    ilist.append(D)
            if len(ilist) > 0:
                Clist.append(sum(ilist)/len(ilist))


    Interlist = []

    for key in clusters:
        for subkey in clusters[key]:
            site1 = subkey

            ilist = []

            for site in dist:

                if site not in clusters[key]:
                    site2 = site
                    D = dist[site1][site2]
                    ilist.append(D)
            if len(ilist) > 0:
                Interlist.append(min(ilist))

    slist = []
    print(len(Clist),len(Interlist))
    for i in range(0,len(Clist)):
        b = Interlist[i]
        a = Clist[i]
        s = (b - a) / max(a,b)
        slist.append(s)
    #print(slist)
    sil = sum(slist)/len(slist)

    return sil

def silhouette_H(clusters,distances):

    #print(clusters)
    #print(len(clusters))
    #for i in clusters:
        #print(clusters[i])
    sil = {}
    for i in range(0,len(clusters)):
        #print(clusters[i])
        dist = distances
        Clist = []
        for key in clusters[i]:
            for subkey in clusters[i][key]:
                site1 = subkey
                ilist = []
                for subkey2 in clusters[i][key]:
                    site2 = subkey2
                    D = dist[site1][site2]
                    ilist.append(D)
                if len(ilist) > 0:
                    Clist.append(sum(ilist)/len(ilist))
        dist2 = distances
        Interlist = []
        for key in clusters[i]:
            for subkey in clusters[i][key]:
                ilist = []
                site1 = subkey
                #print(site1)

                for asd in dist:
                    #print(clusters[i][key])
                    if asd not in clusters[i][key]:
                        site2 = asd
                        #print(site1,site2)
                        D = dist2[site1][site2]
                        ilist.append(D)
                        #print(D)
                if len(ilist) > 0:
                    Interlist.append(min(ilist))

        slist = []
        print(len(Clist),len(Interlist))
        for j in range(0,len(Clist)):
            b = Interlist[j]
            a = Clist[j]
            s = (b - a) / max(a,b)
            slist.append(s)
        sil[i] = sum(slist)/len(slist)
    for kl in sil:
        print(sil)
        print("Iteration {}: {} \n".format(kl,sil[kl]))
    return sil


def centroid(site_a):
    """quick function that calculates the average x, y, and z coordinate for
    an active site's atom coordinates, then subtracts those numbers from the original
    coordinates, thus centering the site around 0."""
    xlist = []
    ylist = []
    zlist = []
    for residue in site_a.residues:
        for atom in residue.atoms:
            xlist.append(atom.coords[0])
            ylist.append(atom.coords[1])
            zlist.append(atom.coords[2])
            #print(atom.coords)
    xmean = sum(xlist) / len(xlist)
    ymean = sum(ylist) / len(ylist)
    zmean = sum(zlist)/ len(zlist)
    xlistnew = []
    ylistnew = []
    zlistnew = []
    for residue in site_a.residues:
        for atom in residue.atoms:
            x_coordn = atom.coords[0] - xmean
            y_coordn = atom.coords[1] - ymean
            z_coordn = atom.coords[2] - zmean
            xlistnew.append(x_coordn)
            ylistnew.append(y_coordn)
            zlistnew.append(z_coordn)
            atom.coords = (x_coordn, y_coordn, z_coordn)


    return site_a


def matrix_create(site_a,site_b):
    """creates an array of all atoms for each residue. """
    V = []
    W = []
    atoms_a = []
    atoms_b = []
    # Align site A and site B
    dict_a,dict_b = sequence_alignment(site_a.residues,site_b.residues)
    res_a_to_add = []
    res_b_to_add = []
    """ Retrieve only the residues that have a match """
    for residue in dict_a:
        if dict_a[residue] != "none":
            res_a_to_add.append(residue)
    for residue in dict_b:
        if dict_b[residue] != "none":
            res_b_to_add.append(residue)
    """ For each residue pair, create a matrix of their atom coordinates.
    In each case, find the one with more less atoms, and only create the matrix
    up to that length, so that the matrices are of equal dimensions. """
    for residue in res_a_to_add:
        index = res_a_to_add.index(residue)
        len1 = len(residue.atoms)
        len2 = len(res_b_to_add[index].atoms)
        for i in range(0,min(len1,len2)):
            atoms_a.append(residue.atoms[i].type)
            x = residue.atoms[i].coords[0]
            y = residue.atoms[i].coords[1]
            z = residue.atoms[i].coords[2]
            V.append(np.asarray([x,y,z],dtype = float))
    for residue in res_b_to_add:
        index = res_b_to_add.index(residue)
        len1 = len(residue.atoms)
        len2 = len(res_a_to_add[index].atoms)
        for i in range(0,min(len1,len2)):
            atoms_b.append(residue.atoms[i].type)
            x = residue.atoms[i].coords[0]
            y = residue.atoms[i].coords[1]
            z = residue.atoms[i].coords[2]
            W.append(np.asarray([x,y,z],dtype = float))
    #print(V)
    #print(W)
    #print(atoms_a)
    #print(atoms_b)
    V = np.asarray(V)
    W = np.asarray(W)
    atoms_a = np.asarray(atoms_a)
    assert(W.shape[0] == V.shape[0])
    return V, W

def find_rotation_matrix(P,Q):
    """ This function finds the rotation matrix U, which rotates an inputed matrix
    to the position which minimizes the RMSD of atom distances. Based on the
    Kabsch algorithm. """

    """ Compute covariance matrix, then use single value decomposition
    function from numpy to get optimal rotation matrix. """
    P,Q = matrix_create(centroid(P),centroid(Q))
    A = np.dot(np.transpose(P),Q)
    V, S, W = np.linalg.svd(A)
    #print("Covariance",A,"V", V,"S", S,"W", W)

    """ If dot product of P and W is negative, correct optimal rotation matrix. """
    if (np.linalg.det(V) * np.linalg.det(W)) < 0:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    U = np.dot(V,W)


    return U,P,Q

def rotate_matrix(P,Q):
    """ Takes rotation matrix U and rotates matrix P relative to matrix Q """
    U,P,Q = find_rotation_matrix(P,Q)

    P = np.dot(P,U)
    #print(P,Q)
    return P,Q




def res_similarity(aa1,aa2):
    """
    Input: two amino acids.
    Output: similarity of the amino acid pair from the similarity matrix.

    The similarity matrix comes from Henikoff, S. and Henikoff, J.G. Amino acid
    substitution matrices from protein blocks. PNAS 1992 89(22):10915.
    The more similar two amino acids are, the better score they will give
    in the alignment. """


    similarity_matrix =[[9,-1,-1,-3,0,-3,-3,-3,-4,-3,-3,-3,-3,-1,-1,-1,-1,-2,-2,-2],
    [-1,4,1,-1,1,0,1,0,0,0,-1,-1,0,-1,-2,-2,-2,-2,-2,-3],
    [-1,1,5,-1,0,-2,0,-1,-1,-1,-2,-1,-1,-1,-1,-1,0,-2,-2,-2],
    [-3,-1,-1,7,-1,-2,-2,-1,-1,-1,-2,-2,-1,-2,-3,-3,-2,-4,-3,-4],
    [0,1,0,-1,4,0,-2,-2,-1,-1,-2,-1,-1,-1,-1,-1,0,-2,-2,-3],
    [-3,0,-2,-2,0,6,0,-1,-2,-2,-2,-2,-2,-3,-4,-4,-3,-3,-3,-2],
    [-3,1,0,-2,-2,0,6,1,0,0,1,0,0,-2,-3,-3,-3,-3,-2,-4],
    [-3,0,-1,-1,-2,-1,1,6,2,0,-1,-2,-1,-3,-3,-4,-3,-3,-3,-4],
    [-4,0,-1,-1,-1,-2,0,2,5,2,0,0,1,-2,-3,-3,-2,-3,-2,-3],
    [-3,0,-1,-1,-1,-2,0,0,2,5,0,1,1,0,-3,-2,-2,-3,-1,-2],
    [-3,-1,-2,-2,-2,-2,1,-1,0,0,8,0,-1,-2,-3,-3,-3,-1,2,-2],
    [-3,-1,-1,-2,-1,-2,0,-2,0,1,0,5,2,-1,-3,-2,-3,-3,-2,-3],
    [-3,0,-1,-1,-1,-2,0,-1,1,1,-1,2,5,-1,-3,-2,-2,-3,-2,-3],
    [-1,-1,-1,-2,-1,-3,-2,-3,-2,0,-2,-1,-1,5,1,2,1,0,-1,-1],
    [-1,-2,-1,-3,-1,-4,-3,-3,-3,-3,-3,-3,-3,1,4,2,3,0,-1,-3],
    [-1,-2,-1,-3,-1,-4,-3,-4,-3,-2,-3,-2,-2,2,2,4,1,0,-1,-2],
    [-1,-2,0,-2,0,-3,-3,-3,-2,-2,-3,-3,-2,1,3,1,4,-1,-1,-3],
    [-2,-2,-2,-4,-2,-3,-3,-3,-3,-3,-1,-3,-3,0,0,0,-1,6,3,1],
    [-2,-2,-2,-3,-2,-3,-2,-3,-2,-1,2,-2,-2,-1,-1,-1,-1,3,7,2],
    [-2,-3,-2,-4,-3,-2,-4,-4,-3,-2,-2,-3,-3,-1,-3,-2,-3,1,2,11]]
    aminoacids = ['CYS','SER','THR','PRO','ALA','GLY','ASN','ASP','GLU','GLN','HIS','ARG','LYS','MET','ILE','LEU','VAL','PHE','TYR','TRP']

    i = aminoacids.index(aa1)
    j = aminoacids.index(aa2)

    s = similarity_matrix[i][j]
    return s




def sequence_alignment(residues_a,residues_b):

    """ This function aligns two amino acid sequences.

    Input: two lists of residue objects (type and number).
    Output: two lists of aligned sequences.

    The goal of using this alignment is to determine which residues
    can be ignored in the rotation matrix and in the calculation of RMSD. """
    reslist_a = [None]
    reslist_b = [None]
    for residue in residues_a:
        reslist_a.append(residue.type)
    for residue in residues_b:
        reslist_b.append(residue.type)
    len1 = len(reslist_a)
    len2 = len(reslist_b)
    F = np.zeros(shape = (len1,len2))
    # define gap penalty "d"
    d = -5
    # Now we create the matrix F. The rows of F represent each residue in one sequence,
    # while the columns represent the other. First, fill in the top row and first column
    # with the gap penalties that would correspond to a non-alignment.
    for i in range(0,len1):
        F[i,0] = d * i
    for j in range(0,len2):
        F[0,j] = d * j
        #Now, for each spot in the matrix, determine whether it would be best to
        # match the two amino acids or have an insert or deletion and place
        # that number in the matrix.
    for i in range (1,len1):
        for j in range(1,len2):
            match = F[i-1,j-1] + res_similarity(reslist_a[i],reslist_b[j])
            delete = F[i-1,j-1] + d
            insert = F[i,j-1] + d
            F[i,j] = max(match,insert,delete)


    alignment_a = []
    alignment_b = []
    dict_a = {}
    dict_b = {}
    i = len1 - 1
    j = len2 - 1

    #Here we walk backwards through the F matrix to determine whether each position
    # on each sequence is optimal with a residue or a gap.
    while i > 0 or j > 0:
        # This walks us backwards diagonally, meaning an amino acid is inserted
        # at the beginning of the alignment list for each sequence. The following
        # should be true if we decided for F[i-1,j-1] to continue to F[i,j]
        # with both amino acids, scoring via similarity function:
        if i > 0 and j > 0 and F[i,j] == F[i-1,j-1] + res_similarity(reslist_a[i],reslist_b[j]):

            alignment_a.insert(0,reslist_a[i])
            alignment_b.insert(0,reslist_b[j])
            dict_a[residues_a[i-1]] = residues_b[j-1]
            dict_b[residues_b[j-1]] = residues_a[i-1]
            i = i - 1
            j = j - 1
        # Inserts a gap for sequence b
        elif i > 0 and F[i,j] == F[i-1,j] + d:
            alignment_a.insert(0,reslist_a[i])
            alignment_b.insert(0,"none")
            dict_a[residues_a[i-1]] = "none"
            #dict_b["none" + str(j - 1)] = residues_a[i-1]
            i = i - 1
        # Inserts a gap for sequence a
        else:
            alignment_a.insert(0,"none")
            alignment_b.insert(0,reslist_b[j])
            #dict_a["none" + str(i-1)] = residues_b[j-1]
            dict_b[residues_b[j-1]] = "none"
            j = j - 1
    """
    print("alignment a",alignment_a)
    print("alignment b",alignment_b)
    print("dict a",dict_a)
    print("dict b",dict_b)
    """
    return (dict_a,dict_b)


def compute_similarity(P,Q):
    P,Q = rotate_matrix(P,Q)
    m = len(P[0])
    n = len(P)
    rmsd = 0
    for i,j in zip(P,Q):
        for k in range(m):
            rmsd = rmsd + sum([(i[k]-j[k])**2])

    rmsd = np.sqrt(rmsd/n)
    return rmsd


def find_distances(active_sites):
    distances = {}
    for site in active_sites:
        site1 = site
        distances[site1] = {}
        for j in active_sites:
            site2 = j
            score = compute_similarity(site1,site2)
            distances[site1][site2] = score
    return distances

def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
            Note: Using K-Medioids (partitioning around medioids).
            The reason I'm using medioids instead of k-means is that
            I don't have a way of calculating the "mean active site" (without
            getting into some "average sequence" reconstruction or something).
            Thus, within each cluster, the goal will be to try out each active site
            as the mean, and see which one works best. This then becomes the new mediod,
            and we re-calculate which cluster every active site belongs to based on
            the new medioids.
    """
    # Fill in your code here!
    #print(random.choice(active_sites))
    distances = find_distances(active_sites)
    n_clusters = int(sys.argv[5])
    num_iterations = int(sys.argv[4])
    min_clustersize = 3
    centers = []
    """choose initial medioid points. The while loop ensures that no duplicates
    are chosen. """
    for i in range (0,n_clusters):
         #make a list of random centers for n_cluster # of centers
        temp = random.choice(active_sites)
        while temp in centers:
            temp = random.choice(active_sites)
        centers.append(temp)
    #print(centers)
    clusters = {}
    """ create dictionary containing n_clusters lists, which will be populated
    with members of that cluster """
    for i in range(0,n_clusters):
        clusters[centers[i]]={}
    #for each active site, compute similarity to each cluster center
    for acs in active_sites:
        site1 = acs
        scorelist = []
        centers = []
        # compare to cluster center from "clusters" dictionary
        """ compare each active site to each central point/medioid """
        for key in clusters:
            site2 = key
            score = distances[site1][site2]
            scorelist.append(score)
            centers.append(key)
            """ This bit shouldn't be necessary, but just in case - if a site is
            being tested against itself, it should go in its own cluster. This
            was a problem before my scoring algorithm was robust. """
            if site1 == site2:
                parent = site2
                break
        #print(scorelist)
        index = np.argmin(scorelist)
        """ assign to cluster - use maximum score, unless the two sites are the same,
        in which case assign active site to its own cluster """
        if site1 != site2:
            parent = centers[index]
        clusters[parent][site1] = score
    medioidlist = []
    """ This series of lists checks whether the current iteration is the same as
    the previous iteration, the i-2 iteration, or the i-3 iteration. """
    check_a = []
    check_b = []
    check_c = []
    check_d = []

    for i in range(0,num_iterations):

        newdict = {}
        for key in clusters:
            """ for each key (medioid), swap to all possible alternative medioids,
            and generate a dictionary of scores for each active site. """
            similarities = clusters[key].values() #gets values from dictionary
            scoresum = {}
            for k in clusters[key]:
                medioid = k
                scorelist = []
                """ for each potential medioid, calculate the similarity of all other
                points in the cluster, make a list, sum the list, and append that sum
                to the dictionary scoresum using the medioid name as the key. """
                for j in clusters[key]:
                    site2 = j
                    score = distances[medioid][site2]
                    scorelist.append(score)
                # compute sum of scores for each potential medioid, and append to
                # dictionary of score sums.
                s = sum(scorelist)
                scoresum[medioid] = s
            # medioid with the best sum of scores is "best" medioid.
            """ this if/else statement adds the best medioid to a new dictionary
            if that subdictionary has more than one element, otherwise it chooses
            a random new active site. This avoids choosing outliers as medioids,
            which can cause some clusters to be a single active site big. Additionally,
            the minimum number of active sites per cluster can be specified here. """
            if len(clusters[key].values()) >= min_clustersize:
                best = min(scoresum.items(), key = lambda x: x[1])[0]
            else:
                best = random.choice(active_sites)
                while best in clusters or best in newdict:
                    best = random.choice(active_sites)
            newdict[best] = {}
            medioidlist.append(best)
            #print("best",best)
        """ This series of lists checks whether the current iteration is the same as
        the previous iteration, the i-2 iteration, or the i-3 iteration. """
        check_a = check_b
        check_b = check_c
        check_c = check_d
        check_d = medioidlist
        if set(check_d) == set(check_c) or set(check_d) == set(check_b) or set(check_d) == set(check_a):
            print("Medioids stable; clustering finished.")
            break
        #print(medioidlist)
        medioidlist = []
        clusters = newdict
        print("Current medioids:",clusters)

        #for each active site, compute similarity to each cluster center
        """ This is the same code as the beginning of the algorithm, but
        since I want to both start and end with this piece of code, I have it inside
        and outside the outermost for loop. """
        for acs in active_sites:
            site1 = acs
            scorelist = []
            centers = []
            # compare to cluster center from "clusters" dictionary
            """ compare to each central point/medioid """
            for key in clusters:
                site2 = key
                score = distances[site1][site2]
                scorelist.append(score)
                centers.append(key)
                if site1 == site2:
                    parent = site2
                    break
            #print(scorelist)
            index = np.argmin(scorelist)
            """ assign to cluster - use maximum score, unless the two sites are the same,
            in which case assign active site to its own cluster """
            if site1 != site2:
                parent = centers[index]
            clusters[parent][site1] = score

    print("Copying pdb files to medioid output folders")
    print("Medioid","Active site")
    """for key in clusters:
        for subkey in clusters[key]:
            print(key,"<--", subkey)
            path = os.path.abspath("data/{}.pdb".format(subkey))
            path2 = os.path.abspath("data/Partitioning/K_=_{}/{}/{}.pdb".format(n_clusters,key,subkey))
            path3 = os.path.abspath("data/Partitioning/K_=_{}/{}".format(n_clusters,key))
            if not os.path.exists(path3):
                os.makedirs(path3)
            copyfile(path, path2)"""

    sil = silhouette_P(clusters,distances)
    print("Silhouette score:", sil)
    return clusters




def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!
    num_iterations = int(sys.argv[4])
    clusters = {}
    distances = find_distances(active_sites)
    ok = float(0)
    clusters[ok] = {}
    j = 0
    for acs in active_sites:
        clusters[ok][j] = [acs]
        j = j + 1

    for i in range(1,len(active_sites)):
        ok = ok + 1
        clusters[ok] = {}
    ok = 0
    biglist = []
    bkeylist = []
    for i in range (0,num_iterations):
        #print(i)
        """ For the first iteration, we need to populate our "clusters" dictionary. """


        #print(clusters[i], "\n")
        """ Here, we are going to compare every cluster to every other cluster.
        This involves iterating through all clusters (for key in clusters[i]),
        then going through the members of that cluster (for val lin clusters[i][key]).
        Nested inside those loops are similar functions for the comparison site.
        For each cluster-cluster comparison, we make a list of all comparisons and
        average their values. Thus, I am using average linkage hierarchical clustering.

        Also note: this function will build several lists.
        keylist: used to prevent scoring of the same two active sites twice.
        sim_means: list of mean similarities between clusters.
        key1: list of keys/clusters, to be indexed in the same way as sim_means.
        key2: list of keys/clusters.
        These last three lists will be used for updating the dictionary of clusters
        during each iteration. """
        keylist = []
        sim_means = []
        key1 = []
        key2 = []
        """ Randomize dictionary key order so that we are not biased towards
        clustering active sites earlier in the list, since more than one
        pair of clusters may have the same average similarity. """
        values1 = list(clusters[ok].keys())
        values2 = list(clusters[ok].keys())

        #print("values1:",values1, "values2:",values2)
        for key in values1:
            keylist.append(key)
            for jkey in values2:
                if jkey not in keylist:
                    for val in clusters[ok][key]:
                        site1 = val
                        sim_inner = []
                        for jval in clusters[ok][jkey]:
                            site2 = jval
                            score = distances[site1][site2]
                            sim_inner.append(score)
                            #print("site1:",site1, "site2:",site2,"score:",
                            #score, "list of similarities for cluster:",sim_inner)
                            #print("cluster1:",key, "cluster2:",jkey)
                        if len(sim_inner) > 0:
                            mean = sum(sim_inner) / len(sim_inner)
                            sim_means.append(mean)
                            #print("mean similarity for cluster:",mean)
                            key1.append(key)
                            key2.append(jkey)
                    #print(sim_means)
        """ Merge the two clusters and delete the redundant one. """
        index = np.argmin(sim_means)

        temp = clusters[ok]

        temp[key1[index]].extend(temp[key2[index]])

        del temp[key2[index]]

        ok = ok + 1
        clusters[ok] = temp
        smalllist = []
        skeylist = []
        for key3 in temp:
            smalllist.append(temp[key3])
            skeylist.append(key3)
        biglist.append(smalllist)
        bkeylist.append(skeylist)
        #print(biglist)
        #print("ok = ",ok)
        #print("clusters ok",clusters[ok])
        #print("temp",temp)

        #print("clusters[%s]:" % i,clusters[i])
        #print("# clusters:", len(clusters[i]))
        if len(clusters[ok]) == 1:
            break

        #print(clusters[ok-1])
    #print(biglist)
    for i in range(0,len(active_sites)):
        r = len(biglist[i])
        clusters[i] = {}
        for j in range(0,r):
            clusters[i][j] = 1
        #print(clusters[ok])
    #sil = silhouette_H(clusters,distances)
    #print("Silhouette scores:",sil)
    return clusters

"""
# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(active_sites)
    write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clusterings = cluster_hierarchically(active_sites)
    print(clusterings)
    write_mult_clusterings(sys.argv[3], clusterings)
"""
"""if sys.argv[1][0:2] == '-S':
    print("Testing similarity metric")"""
""" This script takes every file in 'data' and runs it through compute_similarity with
    every other file, without being redundant (i.e. each pair only goes through
    compute_similarity once) """
"""
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
"""
