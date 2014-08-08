import argparse
import sys
import numpy as np
import libflatclust

# usage: python test_flatclust.py --matrixfile ../../data/rnd_256_256.csv --dictfile ../../data/reuters_dictionary.txt --clusters 16 --infile_W ../../data/flatclust_init_w.csv --infile_H ../../data/flatclust_init_h.csv --miniter 1 --algorithm HALS --maxiter 5000

''' 
if cmp -s "$xml_file" ../../../data/test/flatclust_rnd_clusters_16.xml; then
    echo "XML file test passed"
else
    echo "XML file test failed"
fi
if cmp -s "$assign_file" ../../../data/test/flatclust_rnd_assignments_16.csv; then
    echo "assignment file test passed"
else
    echo "assignment file test failed"
fi
'''

rng = libflatclust.PyRandom()
rng.SeedFromTime()
RNG_CENTER = 0.5
RNG_RADIUS = 0.5

parser = argparse.ArgumentParser()
parser.add_argument("--matrixfile", action="store", required=True,    metavar="matrixfile")
parser.add_argument("--dictfile",   action="store", required=True,    metavar="dictfile")
parser.add_argument("--clusters",   action="store", required=True,    metavar="clusters",   type=int)
parser.add_argument("--algorithm",  action="store", required=False,   metavar="algorithm",    default="BPP")
parser.add_argument("--infile_W",   action="store", required=False,   metavar="infile_W",     default="")
parser.add_argument("--infile_H",   action="store", required=False,   metavar="infile_H",     default="")
parser.add_argument("--tol",        action="store", required=False,   metavar="tol",        type=float,  default=0.0001)
parser.add_argument("--outdir",     action="store", required=False,   metavar="outdir",       default="")
parser.add_argument("--miniter",    action="store", required=False,   metavar="miniter",    type=int,  default=5)
parser.add_argument("--maxiter",    action="store", required=False,   metavar="maxiter",    type=int,  default=5000)
parser.add_argument("--maxterms",   action="store", required=False,   metavar="maxterms",     default=5)
parser.add_argument("--maxthreads", action="store", required=False,   metavar="maxthreads",   default=4)
parser.add_argument("--verbose",    action="store", required=False,   metavar="verbose",      default=True)
parser.add_argument("--format",     action="store", required=False,   metavar="format",       default="XML")
parser.add_argument("--clustfile",  action="store", required=False,   metavar="clustfile",    default="")
parser.add_argument("--assignfile", action="store", required=False,   metavar="assignfile",   default="")
args = parser.parse_args()


opts = libflatclust.PyParseCommandLine2(args)

libflatclust.PyNmfInitialize(len(sys.argv), sys.argv)
print "Initialized: ", libflatclust.PyNmfIsInitialized()
ok = True

# Loading dictionary
dictionary = np.genfromtxt(args.dictfile, dtype="str")
#print dictionary


# Load matrix, either sparse or dense
A = libflatclust.PyDoubleSparseMatrix()
if libflatclust.PyIsSparse(args.matrixfile):
    (ok, m, n, nnz) = libflatclust.PyLoadSparseMatrix(args.matrixfile, A)
    #print ok
elif libflatclust.PyIsDense(args.matrixfile):
    (ok, m, n, buf_a) = libflatclust.PyLoadDenseMatrix(args.matrixfile)
print "Matrix loaded sucessfully", ok

opts["clust_opts"]["nmf_opts"]["height"] = m
opts["clust_opts"]["nmf_opts"]["width"] = n
k = opts["clust_opts"]["nmf_opts"]["k"]
for key in opts.keys():
    print key, opts[key]

ldim_a = m
ldim_w = m
ldim_h = k
height_w = m
width_w = k
height_h = k
width_h = n

# Initializing factor matrices W and H
if not args.infile_W:
    (ok, buf_w) = libflatclust.PyRandomMatrix(m, k, rng, RNG_CENTER, RNG_RADIUS)
else:
    (ok, buf_w) = libflatclust.PyLoadDelimitedFile(height_w, width_w, args.infile_W)
print "Factor matrices W loaded initialized successfully", ok 
if not args.infile_H:
    (ok, buf_h) = libflatclust.PyRandomMatrix(k, n, rng, RNG_CENTER, RNG_RADIUS)
else:
    (ok, buf_h) = libflatclust.PyLoadDelimitedFile(height_h, width_h, args.infile_H)
print "Factor matrices H loaded initialized successfully", ok

stats = libflatclust.PyNmfStats()

if A.Size() > 0:
    (res, buf_w, buf_h) = libflatclust.PyFlatClustSparse( opts, A, buf_w, buf_h, ldim_w, ldim_h, stats)
else:
    (res, buf_w, buf_h) = libflatclust.PyFlatClust(opts, buf_a, buf_w, buf_h, ldim_a, ldim_w, ldim_h, stats)
print "NMF success", res


libflatclust.PyFlatClustAndWriteResults(opts, opts["assignfile"], opts["clustfile"], dictionary, 
                                        buf_h, buf_w, ldim_h, ldim_w, m, n, k, stats)

libflatclust.PyNmfFinalize()