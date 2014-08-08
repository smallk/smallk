# -*- coding: utf-8 -*-
# Copyright 2013,2014 Georgia Institute of Technology.
#
# Licensed under the Apache License, Version 2.0 (the “License”); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by
# applicable law or agreed to in writing, software distributed under the License
# is distributed on as “AS IS” BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the specific language
# governing permissions and limitations under the License.

import argparse
import time
import sys
import numpy as np
import libhierclust

rng = libhierclust.PyRandom()
rng.SeedFromTime()
RNG_CENTER = 0.5
RNG_RADIUS = 0.5

parser = argparse.ArgumentParser()
parser.add_argument("--matrixfile", action="store", required=True,  metavar="matrixfile")
parser.add_argument("--dictfile",   action="store", required=True,  metavar="dictfile")
parser.add_argument("--clusters",   action="store", required=True,  metavar="clusters",   type=int)
parser.add_argument("--infile_W",   action="store", required=False, metavar="infile_W",     default="")
parser.add_argument("--infile_H",   action="store", required=False, metavar="infile_H",     default="")
parser.add_argument("--tol",        action="store", required=False, metavar="tol",        type=float,  default=0.0001)
parser.add_argument("--outdir",     action="store", required=False, metavar="outdir",       default="")
parser.add_argument("--miniter",    action="store", required=False, metavar="miniter",    type=int,  default=5)
parser.add_argument("--maxiter",    action="store", required=False, metavar="maxiter",    type=int,  default=5000)
parser.add_argument("--maxterms",   action="store", required=False, metavar="maxterms",   type=int,  default=5)
parser.add_argument("--maxthreads", action="store", required=False, metavar="maxthreads", type=int,  default=4)
parser.add_argument("--unbalanced", action="store", required=False, metavar="unbalanced", type=float,  default=0.1)
parser.add_argument("--trial_allowance", action="store", required=False, metavar="trial_allowance", default=3)
parser.add_argument("--flat",       action="store", required=False, metavar="flat",         default=0)
parser.add_argument("--verbose",    action="store", required=False, metavar="verbose",      default=True)
parser.add_argument("--format",     action="store", required=False, metavar="format",       default="XML")
parser.add_argument("--treefile",  action="store", required=False, metavar="clustfile",    default="")
parser.add_argument("--assignfile", action="store", required=False, metavar="assignfile",   default="")
args = parser.parse_args()

opts = libhierclust.PyParseCommandLine(args)
libhierclust.PyNmfInitialize(len(sys.argv), sys.argv)
print "Initialized: ", libhierclust.PyNmfIsInitialized()

ok = True

# loading dictionary file
dictionary = np.genfromtxt(args.dictfile, dtype="str")

# loading matrix A, the data matrix
A = libhierclust.PyDoubleSparseMatrix()
if libhierclust.PyIsSparse(args.matrixfile):
    (ok, m, n, nnz) = libhierclust.PyLoadSparseMatrix(args.matrixfile, A)
elif libhierclust.PyIsDense(args.matrixfile):
    (ok, m, n, buf_a) = libhierclust.PyLoadDenseMatrix(args.matrixfile)
print "Matrix loaded sucessfully", ok

num_clusters = opts["clust_opts"]["num_clusters"]
num_initializers = 2*num_clusters
opts["clust_opts"]["nmf_opts"]["height"] = m
opts["clust_opts"]["nmf_opts"]["width"] = n
opts["clust_opts"]["nmf_opts"]["k"] = 2
opts["clust_opts"]["nmf_opts"]["verbose"] = True # comment this line out for silent run
opts["clust_opts"]["nmf_opts"]["max_threads"] = 8
ldim_a = m

# load initializer matrices
height_w = m
width_w = 2
height_h = 2
width_h = n

if not args.infile_W:
    required_size = height_w*width_w
    w_initializers = [[0]*required_size for i in xrange(num_initializers)] 
    for i in xrange(num_initializers):
        (ok, w_initializers[i]) = libhierclust.PyRandomMatrix(height_w, width_w, rng, RNG_CENTER, RNG_RADIUS, required_size)
else:
    (ok, w_initializers) = libhierclust.PyLoadMatrixArray(height_w, width_w, args.infile_W, num_initializers)
print "Factor matrices W loaded initialized successfully", ok 
if not args.infile_H:
    required_size = height_h*width_h
    h_initializers = [[0]*required_size for i in xrange(num_initializers)] 
    for i in xrange(num_initializers):
        (ok, h_initializers[i]) = libhierclust.PyRandomMatrix(height_h, width_h, rng, RNG_CENTER, RNG_RADIUS, required_size)
else:
    (ok, h_initializers) = libhierclust.PyLoadMatrixArray(height_h, width_h, args.infile_H, num_initializers)
print "Factor matrices H loaded initialized successfully", ok 

# run nmf
tree = libhierclust.PyTree()
stats = libhierclust.PyClustStats()
start = time.time()
if A.Size() > 0:
    (ok, buf_w, buf_h, assignments) = libhierclust.PyClustSparse(opts, A, w_initializers, h_initializers, tree, stats, m, n, num_clusters)
else: 
    (ok, buf_w, buf_h, assignments) = libhierclust.PyClust(opts, buf_a, ldim_a, w_initializers, h_initializers, tree, stats, m, n, num_clusters)
print "NMF ran successfully", ok

if opts["clust_opts"]["flat"]:
    k = num_clusters
    assignments_flat = libhierclust.PyComputeAssignments(buf_h, k, k, n)
    term_indices = libhierclust.PyTopTerms(opts["clust_opts"]["maxterms"], buf_w, m, m, k)

end = time.time()
elapsed = end-start
print "Elapsed wall clock time: ", elapsed, "s"

num_converged = stats.nmf_count - stats.max_count
print num_converged, "/", stats.nmf_count, " factorizations converged."

#Outputting files
print "Writing output files...",
ok = libhierclust.PyWriteAssignmentsFile(assignments, opts["assignfile"])
ok = tree.Write(opts["treefile"], opts["format"], dictionary)
if opts["clust_opts"]["flat"]:
    libhierclust.FlatClustWriteResults(opts["outdir"], assignments_flat, dictionary, term_indices, opts["format"], 
                                       opts["clust_opts"]["maxterms"], n, opts["clust_opts"]["num_clusters"])
print "done."

