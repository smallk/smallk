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

import os
import sys
import argparse
import numpy as np
import time

import libpreprocess

infile = "matrix.mtx"
indict = "dictionary.txt"
indocs = "documents.txt"

parser = argparse.ArgumentParser()
parser.add_argument("--indir",         action="store", required=True,  metavar="indir")
parser.add_argument("--outdir",        action="store", required=False, metavar="outdir",        default="./")
parser.add_argument("--docs_per_term", action="store", required=False, metavar="docs_per_term", default=3)
parser.add_argument("--terms_per_doc", action="store", required=False, metavar="terms_per_doc", default=5)
parser.add_argument("--maxiter",       action="store", required=False, metavar="maxiter",       default=1000)
parser.add_argument("--precision",     action="store", required=False, metavar="precision",     default=4)
parser.add_argument("--boolean_mode",  action="store", required=False, metavar="precision",     default=0)

args = parser.parse_args()
# make output path end in "/" if necessary
if args.indir[-1] != "/":
    args.indir += "/"
print "sys.argv = ", sys.argv

pathtofile = args.indir+infile
pathtodict = args.indir+indict
pathtodocs = args.indir+indocs

if args.outdir=="./":
    print "default output directory: current directory"

outfile = args.outdir+"reduced_matrix.mtx"
outdict = args.outdir+"reduced_dictionary.txt"
outdocs = args.outdir+"reduced_documents.txt"

dictionary = np.loadtxt(pathtodict, dtype=str)
documents = np.loadtxt(pathtodocs, dtype=str)
num_docs = len(documents)
num_terms = len(dictionary)

A = libpreprocess.PyDoubleSparseMatrix()
t0 = time.time()
(height, width, nonzeros) = libpreprocess.PyLoadMatrixMarketFile(pathtofile, A)
t1 = time.time()
print "Input file load time: ", t1-t0, "s"

boolean_mode=(args.boolean_mode != 0)
t2 = time.time()
M = libpreprocess.PyDoubleTermFrequencyMatrix(A, boolean_mode)
(term_indices, doc_indices, scores) = libpreprocess.PyPreprocess_tf(M, height, width, args.maxiter, args.docs_per_term, args.terms_per_doc)
t3 = time.time()
print "Processing time: ", t3-t2, "s\n"

print "Writing output matrix " + outfile
t4 = time.time()
libpreprocess.PyWriteMtxFile(M, "./reduced_matrix.mtx", scores, args.precision)
t5 = time.time()
print "Output file write time: ", t5-t4, "s"

print "Writing dictionary file " + outdict
print "Writing documents file " + outdocs
t6 = time.time()
libpreprocess.WriteStringsToFile(outdict, dictionary, term_indices, M.Height())
libpreprocess.WriteStringsToFile(outdocs, documents, doc_indices, M.Width())
t7 = time.time()
print "Dictionary + document write time: ", t7-t6, "s"
