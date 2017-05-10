from argparse import ArgumentParser
from functools import partial
from multiprocessing import Pool
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from Bio import SeqIO
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from dna2vec.dna2vec.multi_k_model import MultiKModel



__path_to_trained = './dna2vec/pretrained/dna2vec-20161219-0153-k3to8-100d-10c-29320Mbp-sliding-Xat.w2v'

##preprocessing related
def read_fasta(filename):
    record = SeqIO.parse(filename, "fasta")
    return record

def read_chromosome(record, grna, mm_tol):
    # check the k of this
    watson = str(record.seq)
    crick = str(record.reverse_complement().seq)

    records = [watson,crick]
    k = len(grna)
    targets=[]
    for record in records:
        record = record.upper()
        chromosome = map(lambda ele: ele[-(k+1):-1], filter(lambda ele: len(ele)>20, record.split('GG')))
        targets += list(filter(lambda x: hamming_distance(grna,x, mm_tol), chromosome))
    return targets

def scan_genome(genome, query, mm_tol,pool = False):
    rc = partial(read_chromosome, grna=query, mm_tol=mm_tol)
    targets =[]
    if pool:
        targets += pool.map(rc,(record for record in genome))
    else:
        targets = np.array(list(map(rc,(record for record in genome))))
    targets= [query]+[item for sublist in targets for item in sublist]
    return targets

def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

def hamming_distance(s1,s2, threshold):
    #compute hamming distance of two sequences
    distance = 0
    for c1,c2 in zip(s1,s2):
        if c1 != c2:
            distance += 1
        if distance > threshold:
            return False
    return distance

def vectorize(targets, model, chunksize = 5):
    ## given targets and a model, compute their vector representation
    vectors = map(lambda x: reversed(list(chunkstring(x,chunksize))), targets)
    output = []
    for lst in vectors:
        sublst = np.array([])
        for ele in lst:
            sublst = np.append(sublst, model.vector(ele))
        sublst = np.reshape(np.array(sublst),(500,))
        output.append(sublst)
    output_vectors = np.array(output)
    return output_vectors 

##kernel functions
def alpha_kernel(distances, a, k):
    num_targets = distances.shape[0]
    sorted_dists = np.argsort(distances)
    i = np.arange(0,num_targets)
    j = sorted_dists[:,k]
    knn_dist = distances[i,j] # bandwidth(x) = distance to k-th neighbor of x
    step = k 
    if any(knn_dist==0):
        print("invalid knn distances encountered due to duplicates. stepping k..")
    while any(knn_dist==0):
        step = int(step*4/3)
        indices = np.where(knn_dist==0)
        print("new k = %d" % (step))
        i = indices
        j = sorted_dists[i,step]
        knn_dist[indices] = distances[i,j]
    pdx = (distances / knn_dist).T 
    g_kernel = np.exp(-1*(pdx**a))
    
    return g_kernel

def dist_func(vectors, mk, metric="minkowski"):
    if metric=="minkowski":
        return squareform(pdist(vectors, metric="minkowski", w=mk))
    else:
        if np.all(mk==0):
            mk[:]=1
        vectors = np.multiply(vectors.T, mk)
        vectors = vectors.T
        return squareform(pdist(vectors, metric=metric))

##diffusion related

def diff_potential(diff_op, t):
    X = np.linalg.matrix_power(diff_op,t) #diffused diffusion operator
    X[X == 0] = np.finfo(float).eps #handling zeros
    X[X <= np.finfo(float).eps] = np.finfo(float).eps #handling small values
    diff_potential = -1*np.log(X) #diffusion potential
    return diff_potential

def diff_op(d):
    #compute diffusion operator from symmetric kernel
    diff_deg = np.diag(np.sum(d,0))
    diff_op = np.dot(np.diag(np.diag(diff_deg)**(-1)),d)
    return diff_deg, diff_op

def affinity_matrix(p):
    #row normalized diffusion affinity matrix used for VNE calculation
    # it can be shown that the eigenvalues of this matrix correspond to the eigenvalues of the diffusion matrix
    rowsums = np.sqrt(np.sum(p,axis=1))
    A = np.ndarray(p.shape)
    for x in range(0,p.shape[0]):
        for y in range(0,p.shape[0]):
            A[x,y] = np.divide(np.multiply(rowsums[x], p[x,y]), rowsums[y])
    A[np.isnan(A)] = 0
    A[np.isinf(A)] = 0
    A[A == 0] = np.finfo(float).eps #handling zeros
    A[A <= np.finfo(float).eps] = np.finfo(float).eps #handling small values
    
    return A

def normalize_eigenvalues(l,t):
    # power eigenvalues and normalize them for VNE computation
    l_t = np.absolute(np.power(l,t))
    l_t[l_t == 0] = np.finfo(float).eps #handling zeros
    l_t[l_t <= np.finfo(float).eps] = np.finfo(float).eps #handling small values
    return l_t/np.sum(l_t)

def VNE(eigen_vals, t):
    #A is the diffusion affinity matrix
    eig_normal = normalize_eigenvalues(eigen_vals, t)
    h = np.multiply(eig_normal, np.log(eig_normal))
    h = -1 * np.sum(h)
    return h

def find_t(p, rng, diff_A = False, sd = 0.1):
    # given a diffusion operator, compute its VNE over various t, then return the first t value with a low second derivative
    if not diff_A:
        A = affinity_matrix(p)
    else:
        A = p

    x, y = np.arange(0,rng+1), np.zeros((rng+1,))
    
    eigensystem = np.linalg.eig(A)
    eig_vals = eigensystem[0]
    for t in x:
        y[t] = VNE(eig_vals,t)

    derivs = approximate_derivative(y,2)
    t = np.where(derivs<sd)[0] #np.where returns a tuple
    # we want the first t below the second derivative gate
    return t[1], np.array([x,y]).T

def approximate_derivative(v,d=1):
    #simple 1-d derivative approximation on a vector
    if d == 1:
        func = lambda x,h,f: (f[x+h] - f[x-h])/(2*h)
    if d == 2:
        func = lambda x,h,f: (f[x+h]-2*f[x]+f[x-h])/(h**2)
    derivs = np.ones(v.shape)
    for i in range(1,len(v)-1):
        derivs[i] = func(i,1,v)
    return derivs

def approximate_derivative(v,d=1):
    #simple 1-d derivative approximation on a vector
    if d == 1:
        func = lambda x,h,f: (f[x+h] - f[x-h])/(2*h)
    if d == 2:
        func = lambda x,h,f: (f[x+h]-2*f[x]+f[x-h])/(h**2)

    derivs = np.ones((len(v),))
    for i in range(1,len(v)-1):
        derivs[i] = func(i,1,v)
    return derivs

def density_function(v):
    return(v/sum(v))

def z_score(v):
    return (v-np.mean(v))/(np.sqrt(np.var(v)))

def normal_pdf(v):
    mu = 0
    sigma_sq = 1
    term1 = 1/np.sqrt(2*np.pi*sigma_sq)
    pdf = lambda x: term1 * np.exp(-1* (((x-mu)**2)/2*sigma_sq))
    return pdf(v)

def exponential_pdf(v):
    l = 1
    pdf = lambda x: l*np.exp((-1/l)*x)

    return pdf(v)

def feature_scale(v):
    mi = np.min(v,0)
    ma = np.max(v,0)
    return ((v-mi)/(ma-mi))
if __name__ == "__main__":
    parser = ArgumentParser(description="MAGE- Markov Affinity Guide rna Extraction")
    parser.add_argument('-sgr', dest='sgrna', required=True, type=str, nargs=1, help="FASTA of sgrna")
    parser.add_argument('-gnm', dest='genome', required=True, type=str, nargs=1, help="Genome to target")
    parser.add_argument('-hd', dest='hd', required=False, type=int, nargs=1, default=[5], help="Hamming distance to tolerate across offtargets")
    parser.add_argument('-mk', dest='mk', required=False, type=int, nargs=1, default=[15], help="Integer weighting to apply to the seed region")
    parser.add_argument('-k', dest='k', required=False, type=int, nargs=1, default=[50], help="k-nn adaptive bandwidth for kernel")
    parser.add_argument('-a', dest='a', required=False, type=int, nargs=1, default=[5], help="Falloff exponent of alpha decaying kernel")
    parser.add_argument('-sd', dest ='sd', required=False, type=float, nargs=1, default=[0.1], help="Second derivative target for VNE")
    parser.add_argument('-t', dest='t', required=False, type=int, nargs=1, default=[False], help="Powering to generate diffusion operator.  Default (0) uses 2nd derivative based optimization")
    parser.add_argument('-trange', dest='trange', required=False, type=int, nargs=1, default=[20], help="range of t values to optimize on, default is 20")
    parser.add_argument('-out', dest='output', required=True, type=str, nargs=1, help='Output directory')
    parser.add_argument('-mds', dest='mds', required = False, type=int, nargs = 1, default=[0], help='compute & plot mds embedding, 0: No mds; 1: classical mds; 2: nonmetric mds')
    parser.add_argument('-vne', dest='vne', required=False, type=int, nargs=1, default=[0], help='save & plot von neumann entropy')
    parser.add_argument('-num', dest='num', required=False, type=int, nargs=1 , default=[0], help='number of top hits to retain, affects final pdists; 0: retains all')
    parser.add_argument('-cores', dest='cores', required=False, type=int, nargs=1, default=[1], help='number of cores to use, for scanning genome and MDS embeddings.')
    args = parser.parse_args()

    #system params
    cores = args.cores[0]
    #kernel/distance params
    __a = args.a[0]
    __k = args.k[0]
    __hd = args.hd[0]
    __mk = args.mk[0]
    

    #vne params
    __vneflag = args.vne[0]
    __t = args.t[0]
    __trange = args.trange[0]
    __sd = args.sd[0]

    #inputs
    sgr_f = args.sgrna[0]
    gnm_f = args.genome[0]
    #output directory
    output_dir = args.output[0]
    #targets to keep
    final_num = args.num[0]
    #plotting parameters
    __mdsflag = args.mds[0]

    
    print("starting model")
    vec_model = MultiKModel(__path_to_trained)
    print("model started, importing sequences")
    grna_seq = next(read_fasta(sgr_f)).seq
    genome = read_fasta(gnm_f)
    print("scanning genome")
    if cores>1:
        p = Pool(cores)
        targets = scan_genome(genome,str(grna_seq), __hd, p)
        p.close()
    else:
        targets = scan_genome(genome,str(grna_seq), __hd)

    print("vectorizing genome")
    vecs = vectorize(targets, vec_model, chunksize=4)
    #set up weighting
    mk = np.ones(vecs.shape[1])
    n_kmers = (vecs.shape[1]/100)
    mk[:(3/n_kmers)*100] = __mk

    #compute distances
    dists = dist_func(vecs,mk)
    print("computing kernel")

    if __k*1.5 > len(targets):
        __k = int(1/10*(len(targets)))
        print("k is too close to the number of targets, tuning k to %d" %__k)
    try:
        kernel = alpha_kernel(dists, __a, __k)
    except:
        print("attempted alpha kernel with %d a and %d k, repeating with known stable values" % (__a, __k))
        kernel = alpha_kernel(dists, 3, len(targets)//2)

    sym_kern = kernel+kernel.T
    print("diffusion operator")
    _, operator = diff_op(sym_kern)
    print("computing VNE")
    if not __t:
        t, h_t = find_t(operator, __trange, __sd)
    else:
        t = __t

    print("finding potentials")
    potentials = diff_potential(operator, t)
    print("potential distances obtained")
    potential_distances = squareform(pdist(potentials))

    scaled_distances = feature_scale(potential_distances)
    d_from_guide = scaled_distances[:,0]
    unscaled_d_from_guide = potential_distances[:,0]
    #output related functions
    if __mdsflag != 0:
        embedder = MDS(metric=(__mdsflag!=2), n_components=2,n_jobs=cores, dissimilarity = 'precomputed')
        mdsembedding = embedder.fit_transform(scaled_distances+scaled_distances.T)
        np.savetxt(output_dir+'/potential_embedding.txt', mdsembedding)
        mdsfig = plt.figure()
        plt.scatter(mdsembedding[1:,0], mdsembedding[1:,1],color='blue')
        plt.scatter(mdsembedding[0,0], mdsembedding[0,1], color='orange')
        mdsfig.savefig(output_dir+"/mds_plot.pdf")
        mdsfig=0

    if __vneflag != 0:
        vnefig = plt.figure()
        plt.plot(h_t[:,1])
        vnefig.savefig(output_dir+"/vne_plot.pdf")
        vnefig=0

    np.savetxt(output_dir+"/potential_distances.txt", potential_distances)
    np.savetxt(output_dir+"/scaled_distance.txt", scaled_distances)
    np.savetxt(output_dir+"/vne_pts.txt", h_t)

    if final_num == 0:
        final_num = len(d_from_guide)
    f_d = d_from_guide[np.argsort(d_from_guide) < final_num]
    u_f_d = unscaled_d_from_guide[np.argsort(d_from_guide) < final_num]

    normal = normal_pdf(f_d)
    expo = exponential_pdf(f_d)
    z = z_score(f_d)
    normalized = density_function(f_d)

    out_array = np.array([np.array(targets),f_d, density_function(normal), density_function(expo), z,u_f_d]).T
    out_df = pd.DataFrame(out_array, columns=["seq", "scaled diffusion distance", "normal distribution probability", "exponential distribution","z score", 'actual diffusion distance'])
    sorted_df = out_df[:1].append(out_df[1:].sort('scaled diffusion distance'))
    sorted_df.to_csv(output_dir+"/offtargets.csv")

