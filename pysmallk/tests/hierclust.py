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

#-----------------------------------------------------------------------------
try:
	import sys
	import numpy as np #must be imported prior to use of pysmallk library
	import pysmallk
except ImportError:
	print 'ImportError: hierclust test failed'
	raise

h = pysmallk.Hierclust()

# This program can be used as a command line tool to run hierarchical NMF clustering.
# It also demonstrates how to use each of the hierclust functions, which can be 
# integrated into other code.

# use the parser function to parse the command line arguments from the user
args = h.parser()

# load in the user-specified matrix and dictionary
h.load_matrix(filepath=args.matrixfile)
h.load_dictionary(dictfile=args.dictfile)

# use the user-provided inputs to run hierarchical clustering
# all keyword arguments are optional 


h.cluster(args.clusters, initdir=args.initdir,
	min_iter=args.miniter, max_iter=args.maxiter, tol=args.tol,
	verbose=args.verbose, trial_allowance=args.trial_allowance,
	unbalanced=args.unbalanced, flat=args.flat, maxterms=args.maxterms,
	max_threads=args.maxthreads)

# write files to system using user-provided filenames
# assignfile and treefile should not include file extensions, e.g. 'tree_output'
# the library will append the appropriate ending to the filename based on the 'format' parameter
h.write_output(args.assignfile, args.treefile, args.fuzzyfile, outdir=args.outdir, format=args.format)

# always call finalize() when done with analysis; this helps clean up system variables
h.finalize()
