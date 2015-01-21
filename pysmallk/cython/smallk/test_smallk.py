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

# example: python test_smallk.py --k 5 --algorithm HALS --data_dir data

import sys
import argparse
import numpy # must import numpy before importing the .so files
import libsmallkpy as smallklib


parser = argparse.ArgumentParser(description="Run NMF via python binding")
parser.add_argument('--matrixfile', action='store', required=True, metavar='matrixfile')
parser.add_argument('--k', action='store', required=True, metavar='kval')
parser.add_argument('--algorithm', action='store', required=False, default='BPP', metavar='algorithm', choices=['MU','RANK2','HALS','BPP'])
parser.add_argument('--tol', action='store', required=False, metavar='tol_val', default=0.005)
parser.add_argument('--infile_W', action='store', required=False, metavar='infile_W', default="")
parser.add_argument('--infile_H', action='store', required=False, metavar='infile_H', default="")
parser.add_argument('--outdir', action='store', required=False, metavar='outdir', default=".")
parser.add_argument('--outprecision', action='store', required=False, metavar="outprecision", default=6)
parser.add_argument('--maxiter', action='store', required=False, metavar="maxiter", default=5000)
parser.add_argument('--miniter', action='store', required=False, metavar="miniter", default=5)
parser.add_argument('--maxthreads', action='store', required=False, metavar="maxthreads", default=8)

# The following parser options are not exposed, since their inputs do not effect change in the program execution
#parser.add_argument('--tolcount', action='store', required=False, metavar='tolcount', default=1)
#parser.add_argument('--stopping', action='store', required=False, metavar='stopping', default='PG_RATIO', choices=['PG_RATIO','DELTA'])
#parser.add_argument('--verbose', action='store_true', required=False, default=False, metavar="verbose", default=1)
#parser.add_argument('--normalize', action='store_true', required=False, default=False, metavar="normalize", default=1)

args = parser.parse_args()

# Initialize NMF
smallklib.py_initialize(len(sys.argv),sys.argv)

# Check Initialization
if (smallklib.py_isInitialized() == False):
    print "Error in Initialization of NMFlib"
    exit(-1)
else: 
    print "Initialization oK"

algorithm = smallklib.get_algorithm(args.algorithm);
k = int(args.k)

#include user inputs, or defaults
smallklib.py_loadMatrix(args.matrixfile)
smallklib.PySetMinIter(args.miniter)
smallklib.PySetOutputPrecision(args.outprecision)
smallklib.PySetMaxIter(args.maxiter)
smallklib.PySetMaxThreads(args.maxthreads)
smallklib.PySetNmfTolerance(args.tol)
smallklib.PySetOutputDir(args.outdir)

smallklib.py_nmf(k, algorithm, args.infile_W, args.infile_H)

smallklib.py_finalize()
